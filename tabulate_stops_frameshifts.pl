#!/usr/bin/perl -w
use strict;
use Bio::Perl;
use Bio::SeqIO;

print STDERR "perl ~/bin/tabulate_stops_frameshifts.pl <protein fasta> <gff> <exonerate output>\n";

my ($fasta, $gff, $exonerate) = @ARGV;
#my ($fasta, $gff, $exonerate) = ("/scratch/ykelkar2/muni/EVM/MELL.evm.out.pepp.cds.fasta","/scratch/ykelkar2/rap/annotation/EllusmuniEVM-rapGap-exonerate2.outs.gff3","rapEVM-muniGap-exonerate2.stops_frameshifts.txt");
#my ($fasta, $gff, $exonerate) = ("/scratch/ykelkar2/muni/EVM/MELL.evm.out.pepp.cds.fasta","/scratch/ykelkar2/muni/annotation/muni_cufflinks_male/EllusrapEVM-muniGap-exonerate2.outs.gff3","ellusEVM-muniGap-exonerate2.stops_frameshifts.txt");

#my ($fasta, $gff, $exonerate) = ("/scratch/ykelkar2/muni/EVM/MRAP.evm.out.pepp.cds.fasta",
#"/scratch/ykelkar2/muni/annotation/muni_cufflinks_male/EllusrapEVM-muniGap-exonerate2.outs.gff3",
#"/scratch/ykelkar2/muni/EVM/cross_annotations/rapEVM-muniGap-exonerate2.stops_frameshifts.txt");
#my ($fasta, $gff, $exonerate) = ("/scratch/ykelkar2/muni/EVM/MRAP.evm.out.pepp.cds.fasta","/scratch/ykelkar2/ellus/annotation/rapmuniEVM-ellusGap-exonerate2.outs.gff3","/scratch/ykelkar2/muni/EVM/cross_annotations/rapEVM-ellusGap-exonerate2.stops_frameshifts.txt");


my $in  = Bio::SeqIO->new(-file => "$fasta" , '-format' => 'Fasta');
open ("GFF", "<$gff") or die;
open ("IN", "<$exonerate");

my %qentries=();
my %tentries=();

while ( my $seq_obj = $in->next_seq ) {
        my $spec = $seq_obj->id;
        my $len = length($seq_obj->seq());
	my $seq = $seq_obj->seq();
	$qentries{$spec}{'len'}=$len/3;
	my $n = 3;
        my @codons = unpack "a$n" x (length( $seq ) /$n ), $seq;
	my $count=0;
	foreach my $codon (@codons){
		$count++;
		$qentries{$spec}{'seq'}{$count}=$codon;
	}
#	print STDERR "$spec, \n";
}

print STDERR "Read codons\n"; 

while (my $line = <GFF>){
	my ($c, $e, $type,$start,$end,$dot1,$strand,$dot2,$attr) = split(/\t/,$line);
	next if $type ne "gene";
#	print STDERR $line;
	$attr=~ /ID=([A-Za-z0-9\-\._]+)/;
	my $id = $1;
	my $source  = $1; $source =~ s/^[0-9]+__//g;
#	print STDERR "$id, source = $source\n";

	next if !exists $qentries{$source};
	$tentries{$c}{join("\t",$start, $end)} = $id;
#	print STDERR "$c, $start, $end,  $id, $source\n";

}

print STDERR "Read GFF\n";
#<STDIN>;

my $q = ();
my $t = ();
my $qstart = ();
my $qend = ();
my $qcov = ();
my $revcomp=0;
my $tstart = ();
my $tend = ();
	
my $importname = ();
	
my @insertstops_starts = ();
my @insertstops_stopcodon = ();

my @replacestops_starts = ();
my @replacestops_stopcodon = ();
my @replacestops_AminoAcid = ();
my @replacestops_qcodon = ();

my @frameshifts_starts = ();
my @frameshifts_lens = ();
my @frameshifts_inserts = ();
my @sequence = ();
my $ncount=0;
while (my $line = <IN>){
	print STDERR $line;
	chomp $line;
	$line =~ /^\s+/;
	$ncount++;

	if ($line =~ /Query:/){
		
#		if (scalar(@frameshifts_starts > 0 or @replacestops_starts > 0 or @insertstops_starts > 0)){
		my $joint= join("\t",join(",",@insertstops_starts),join(",",@insertstops_stopcodon),join(",",@replacestops_starts),join(",",@replacestops_stopcodon),join(",",@replacestops_AminoAcid),join(",",@replacestops_qcodon),join(",",@frameshifts_starts),join(",",@frameshifts_lens),join(",",@frameshifts_inserts),join(",",@sequence)  );
		my $mark = 0;
		$mark = 1 if $joint =~ /[A-Za-z0-9]/;
		print join("\t", $q, $t, $qentries{$q}{'len'}, $qstart, $qend, $tstart, $tend, $importname, $mark, join(",",@insertstops_starts),join(",",@insertstops_stopcodon),join(",",@replacestops_starts),join(",",@replacestops_stopcodon),join(",",@replacestops_AminoAcid),join(",",@replacestops_qcodon),join(",",@frameshifts_starts),join(",",@frameshifts_lens),join(",",@frameshifts_inserts),join(",",@sequence)  ),"\n" if $ncount > 0;
			#print STDERR join("\t", $q, $t, $qstart, $qend, $tstart, $tend, $importname, join(",",@insertstops_starts),join(",",@insertstops_stopcodon),join(",",@replacestops_starts),join(",",@replacestops_stopcodon),join(",",@replacestops_AminoAcid),join(",",@replacestops_qcodon),join(",",@frameshifts_starts),join(",",@frameshifts_lens),join(",",@frameshifts_inserts) ,join(",",@sequence) ),"\n";		
			#<STDIN>;
#		}
		
		$q = ();
		$t = ();
		$qstart = ();
		$qend = ();
		$qcov = ();
		$revcomp=0;
		$tstart = ();
		$tend = ();
		$importname = ();
		@insertstops_starts = ();
		@insertstops_stopcodon = ();
		@replacestops_starts = ();
		@replacestops_stopcodon = ();
		@replacestops_AminoAcid = ();
		@replacestops_qcodon = ();
		@frameshifts_starts = ();
		@frameshifts_lens = ();
		@frameshifts_inserts = ();	
		@sequence = ();
		$line =~ /:\s*([\S]+)/;
		$q = $1;	
	}
	elsif ($line =~ /Query range:/){
		$line =~ /:\s*([0-9]+)\s+\->\s+([0-9]+)/;
		$qstart = $1;		
		$qend = $2;
		$qcov = $qend-$qstart+1;
#		print STDERR "$q $qstart $qend $qcov\n"; #<STDIN>;
	}
	elsif ($line =~ /Target:/){
		$line =~ /:\s*([\S]+)/;
		$t = $1;	
		$revcomp	 = 1 if $line =~ /revcomp/;
	}
	elsif ($line =~ /Target range:/){
		$line =~ /:\s*([0-9]+)\s+\->\s+([0-9]+)/;
		if ($revcomp == 1){$tstart = $2; $tend = $1}
		elsif ($revcomp == 0){$tstart = $1; $tend = $2}
		my @importcoords= (keys  %{$tentries{$t}});
		foreach my $coords (@importcoords){
#			print STDERR $coords, "\n";
			my ($istart, $iend) = split(/\t/,$coords);
			$importname=$tentries{$t}{$coords} if $istart = $tstart+1 and $iend == $tend;
		}
#		print STDERR "$t $tstart $tend... $importname\n"; <STDIN> if $importname !~ /[a-zA-Z0-9]/;
	}
	else{
		my ($qstate, $tstate, $qpos_start,$qpos_end, $tpos_start,$tpos_end, $info1, $info2) = split(/\t/,$line);
		next if $qstate eq "*" and $qpos_end == $qend;
		if ($info1 =~ /frameshift open/){
			push @frameshifts_starts, $qpos_start;
			push @frameshifts_lens, length($tstate);
			push @frameshifts_inserts, $tstate;
			push @sequence, "F";
#	print STDERR "@frameshifts_starts, @frameshifts_lens, @frameshifts_inserts, @insertstops_starts, @insertstops_stopcodon, ",
#	"@replacestops_starts, @replacestops_AminoAcid, @replacestops_stopcodon, @replacestops_qcodon\n";
	#<STDIN>;
		}
		elsif ($info1 =~ /frameshift close/){
			$frameshifts_inserts[$#frameshifts_inserts].=$tstate if $tstate =~ /[A-Za-z]/;
			next;
		}
		if ($qstate eq "*" and ($tstate eq "TAA" or $tstate eq "TGA" or $tstate eq "TAG" )   ){
			push @insertstops_starts , $qpos_start;
			push @insertstops_stopcodon, $qstate;
			push @sequence, "I";
#	print STDERR "@frameshifts_starts, @frameshifts_lens, @frameshifts_inserts, @insertstops_starts, @insertstops_stopcodon, ",
#	"@replacestops_starts, @replacestops_AminoAcid, @replacestops_stopcodon, @replacestops_qcodon\n";
	#<STDIN>;
		}
		if ($qstate ne "*" and ($tstate eq "TAA" or $tstate eq "TGA" or $tstate eq "TAG" )   ){
			push @replacestops_starts, $qpos_start;
			push @replacestops_AminoAcid, $qstate;
			push @replacestops_stopcodon, $tstate;
			push @replacestops_qcodon  , $qentries{$q}{'seq'}{$qpos_start};
			push @sequence, "R";
			
#	print STDERR "@frameshifts_starts, @frameshifts_lens, @frameshifts_inserts, @insertstops_starts, @insertstops_stopcodon, ",
#	"@replacestops_starts, @replacestops_AminoAcid, @replacestops_stopcodon, @replacestops_qcodon\n";
	#<STDIN>;
		}		
	}

}
