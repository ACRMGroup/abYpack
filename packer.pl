#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    abYmod
#   File:       packer.pl
#   
#   Version:    V1.0
#   Date:       14.03.22
#   Function:   Pack chains together at required packing angle based on
#               the best template
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

UsageDie() if(defined($::h) || (scalar(@ARGV) != 1) ||
              (!defined($::angle) && !defined($::a)));

$::angle = $::a if(defined($::a));

# Preparation
my $angle         = $::angle;
my $pdbFile       = shift(@ARGV);
my $referenceFile = '';
my $tmpDir        = util::CreateTempDir("abYpack");
util::Die("Cannot create temp directory") if($tmpDir eq '');
$::q = 1;

# Actual run
my $model = RunPacker($tmpDir, $pdbFile, $angle);
system("cat $model");

# Cleanup
unlink($tmpDir);

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
sub RunPacker
{
    my($tmpDir, $pdbFile, $angle) = @_;
    my $bestError   = 1000.0;
    my $bestModel   = '';
    my $modelNumber = 1;

    # Set a range of angles in the library files. This exploits the
    # roughly normal distribution to use smaller ranges when we are
    # close to the peak of the distribution
    my $range = 0.2;
    if(($angle < -55) || ($angle > -35))
    {
        $range = 7;
    }
    elsif(($angle < -57) || ($angle > -37))
    {
        $range = 5;
    }
    elsif(($angle < -50) || ($angle > -40))
    {
        $range = 2;
    }
    elsif(($angle <= -49) || ($angle > -42))
    {
        $range = 1;
    }
    elsif(($angle < -48) || ($angle > -44))
    {
        $range = 0.5;
    }

    my($light, $heavy) = SplitPDB($tmpDir, $pdbFile);
    if(($light eq '') || ($heavy eq ''))
    {
        util::Die("Could not split PDB file into light and heavy chains\n");
    }

    if(open(my $fh, '<', $config::angleList))
    {
        while (my $row = <$fh>) 
        {
            chomp $row;
            $row =~ s/\#.*//;   # Remove comments
            $row =~ s/^\s+//;   # Remove leading whitespace
            if(length($row))
            {
                my($referenceAngle, $referenceFile) = split(/\s+/, $row);
                if(abs($referenceAngle - $angle) < $range)
                {
                    my $model = FitChainsToReference($tmpDir, $light, $heavy,
                                                     $referenceFile, $modelNumber);
                    my $modelAngle = CalcModelAngle($model);
                    my $error = abs($modelAngle - $angle);
                    if($error < $bestError)
                    {
                        $bestError = $error;
                        $bestModel = $model;
                    }
                    $modelNumber++;
                }
            }
        }
        close $fh;
    }
    else
    {
        util::Die("Could not open reference file with packing angles and filenames: '$config::angleList'\n");
    }

    return($bestModel);
}

#*************************************************************************
#my $modelAngle = CalcModelAngle($model);
sub CalcModelAngle
{
    my($model) = @_;
    my $exe = "abpackingangle -q $model";
    my $angle = util::RunCommand($exe);
}


#*************************************************************************
# my($light, $heavy) = SplitPDB($tmpDir, $pdbFile))
sub SplitPDB
{
    my($tmpDir, $pdbFile) = @_;
    my $light   = GetChain($pdbFile, 'L', $tmpDir);
    my $heavy   = GetChain($pdbFile, 'H', $tmpDir);
    return($light, $heavy);
}

#*************************************************************************
#> $model = FitChainsToReference($tmpDir, $light, $heavy, $referenceFile,
#                                $modelNumber)
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
#  14.03.22   Original based on abYmod  By: ACRM
sub FitChainsToReference
{
    my($tmpDir, $light, $heavy, $referenceFile, $modelNumber) = @_;
    
    my $outfile = "$tmpDir/model_$modelNumber.pdb";

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

packer V1.0 (c) 2020, UCL, Prof. Andrew C.R. Martin

Usage: ./packer.pl [-v=nnn] -angle=angle in.pdb > out.pdb
       -v verbose mode
       -angle|-a  - specify the angle REQUIRED!

Packs an antibody VH/VL to (approximately) the required packing angle
by fitting light and heavy chains to a template    

__EOF

   exit 0;
}

