#!/usr/bin/perl

# Script written by Greg Pearce, Jul 2002
# Usage: geom2xyz file.geom    .OR.    geom2xyz file.md 
# Output: file.xyz
# Converts .geom or .md files from CASTEP into .xyz for animation in xmol 

$bohr = 0.5291772083; 
$filename = $ARGV[0];

open(GEOM, "$filename") || die ("File '$filename' not found!\n");
@geom_content=<GEOM>;
close(GEOM);

$filename =~ s/(\.geom|\.md)/\.xyz/i;  

$answer="";
if (-e($filename)) { 
    while ($answer !~ m/^(y|n)/i ) {
	print STDOUT "Replace existing $filename [yn]?";
	$answer=<STDIN>;
	if ($answer !~ m/^y/i) {
            print STDOUT "Enter new file name :";
            $filename=<STDIN>;
            if (-e($filename)) {
	       die ("File $filename exists.  Exiting.");
            }
	}
    }
}

open(XYZ, ">$filename");

$n=0;
foreach $line (@geom_content) {
    $line =~ /.{78}(.{1})/;
    if ($1 =~ m/^R$/) {
	$n++;
    }
    elsif ($n != 0) {last}
}

$newentry = 1;
$past_header = 0;
foreach $line (@geom_content) {
    $line =~ /.(.{2}).{11}(.{16}).{5}(.{16}).{5}(.{16}).{6}(.{1})/;   
    $atom = $1;
    $x = $2*$bohr;
    $y = $3*$bohr;
    $z = $4*$bohr;
    $info = $5;
    if ($line =~ m/\d+/) {$past_header=1}
    if ($info =~ m/^R$/) {
	if ($newentry == 1) {
	    print XYZ $n, "\n\n";
	    $newentry = 0;
	}
       print XYZ $atom, " ", $x, " ", $y, " ", $z, "\n";
    }
    elsif (($line !~ m/\w+/) && ($past_header == 1)) {$newentry=1}
}









