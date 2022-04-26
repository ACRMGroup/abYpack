#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    abYpack
#   File:       abypack.pl
#   
#   Version:    V1.0
#   Date:       14.03.22
#   Function:   Build an antibody model from a specified set of templates
#   
#   Copyright:  (c) UCL, Prof. Andrew C. R. Martin, 2022
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#   
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0   22.03.22  Extracted from abYmod
#
#*************************************************************************
use strict;

# Add the path of the executable to the library path
use FindBin;
use lib $FindBin::Bin;
use config;
use util;
#use abymod;

UsageDie() if(defined($::h) || (scalar(@ARGV) != 1) || (!defined($::angle) && !defined($::a)));

$::angle = $::a if(defined($::a));

# Preparation
my $angle         = $::angle;
my $pdbFile       = shift(@ARGV);
my $referenceFile = '';
my $tmpDir        = util::CreateTempDir("abYpack");
util::Die("Cannot create temp directory") if($tmpDir eq '');

# Actual run
$referenceFile = FindReferenceFile($angle);
print STDERR "Reference file: $referenceFile\n\n"   if($::v >= 1);
print STDERR "Fitting Light and Heavy chains to reference..." if($::v >= 1);
my $model = FitChainsToReference($tmpDir, $pdbFile, $referenceFile);
print STDERR "done\n" if($::v >= 1);
system("cat $model");

# Cleanup
#unlink($tmpDir);

#*************************************************************************
#> $referenceFile = FindReferenceFile($angle)
#  ------------------------------------------
#  Inputs:   real   $angle      A predicted packing angle
#  
#  Reads the reference angle list and finds the antibody file that best 
#  matches the predicted VH/VL packing angle.
#
#  29.11.16  Original   By: ZCL
#  30.11.16  Tidied up  By: ACRM
sub FindReferenceFile
{
    my($angle) = @_;
    my $referenceFile  = '';
    my $referenceAngle = -9999;
    my $foundFile      = 0;

    if(open(my $fh, '<', $config::angleList))
    {
        while (my $row = <$fh>) 
        {
            chomp $row;
            $row =~ s/\#.*//;   # Remove comments
            $row =~ s/^\s+//;   # Remove leading whitespace
            if(length($row))
            {
                ($referenceAngle, $referenceFile) = split(/\s+/, $row);
                if($angle <= $referenceAngle) 
                {
                    $foundFile = 1;
                    last;
                }
            }
        }
        close $fh;
    }
    else
    {
        util::Die("Could not open reference file with packing angles and filenames: '$config::angleList'\n");
    }

    $referenceFile = '' if(!$foundFile);
    return ($referenceFile);
}

#*************************************************************************
#> $model = FitChainsToReference($tmpDir, $pdbFile, $referenceFile)
#  -------------------------------------------------------------
#  Inputs:   string   $tmpDir          Temporary working directory 
#                                      (already created)
#            \hash    $hBestTemplates  Reference to hash containing the
#                                      selected templates for the two 
#                                      chains. Keyed by chain (L or H)
#            string   $referenceFile   Reference file to which to fit
#  Returns:  string                    Full path to resulting model
#
#  Fits the light and heavy chains together using external program.
#
#  14.0322   Original based on abYmod  By: ACRM
sub FitChainsToReference
{
    my($tmpDir, $pdbFile, $referenceFile) = @_;
    
    my $light   = GetChain($pdbFile, 'L', $tmpDir);
    my $heavy   = GetChain($pdbFile, 'H', $tmpDir);
    my $outfile = "$tmpDir/model.pdb";

    $referenceFile = "$config::abpdblib/$referenceFile";

    my $exe = "$config::bindir/fitlhpdb -s";
    $exe .= " -t $referenceFile";
    
    $exe .= " L:$light";
    $exe .= " H:$heavy";
    $exe .= " $outfile";
    
    util::RunCommand($exe);

    return($outfile);
}


#*************************************************************************
sub GetChain
{
    my($pdbFile, $chain, $tmpDir) = @_;
    my $outFile = "$tmpDir/${chain}.pdb";
    my $exe = "$config::bindir/pdbgetchain $chain $pdbFile $outFile";
    util::RunCommand($exe);
    util::CheckAndDie($outFile, 0, "Cannot create temp file: $outFile");
    return($outFile);
}

#*************************************************************************
#> void UsageDie()
#  ---------------
#  Prints a usage message and exits
#
#  19.09.13  Original  By: ACRM
sub UsageDie
{
    print STDERR <<__EOF;

abYpack V1.0 (c) 2020, UCL, Prof. Andrew C.R. Martin

Usage: ./abypack.pl angle in.pdb > out.pdb

__EOF

   exit 0;
}

