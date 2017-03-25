# Search_pseudogenes

## What this tool does:

My tool tries to identify multiple causes of pseudogenization, and their information is contained in different columns (column names italicized):

A. Insertion of stop codons
- *StopCodon_insertion_starts*
- *StopCodon_insertion_stops*
these are coordinates of insertion sites.

B. Replacement of an amino acid codon with a stop codon
- *StopCodon_replace_starts* -- coordinate of replacement's start position
- *StopCodon_replacement* -- sequence of the new stop codon
- *Original_AminoAcid*    -- identity of the original amino acid
- *Original_codon* -- sequence of the original codon.

C. Frameshift
- *Frameshift_start*    --
- *Frameshift_length*    --
- *Frameshift_insertion*  -- all these are intuitive

The last column (named *Sequence_of_events*) is list the sequence of events, from the starting codon to the last codon of the alignment. As you can see, in many of your cases, each of the pseudogene-like hits has multiple issues.. usually starting with a replacement by stop codon (R) that is often followed by frameshifts (F)

The column *modified* is 1 if above columns are populated, and 0 if not.

The column *t1* is useless. 



## Getting inputs ready

Before running this tool, you need to run exonerate, using protein sequences as queries and a genome as a target. Use the following exonerate parameters:

```
exonerate --percent 60 --model protein2genome --showquerygff yes --showtargetgff yes -q QUERY_PROTEINS.fasta -t TARGET_GENOME.fasta --ryo "RYO\t%qi\t%ti\t%ql\t%tl\t%qal\t%qab\t%qae\t%tal\t%tab\t%tae\t%et\t%ei\t%es\t%em\t%pi\t%ps\t%g\nTransitionStart\n%V{%Pqs\t%Pts\t%Pqb\t%Pqe\t%Ptb\t%Pte\t%Pn\t%Pl\n}TransitionEnd\nTargetSeq\n%qs\nAligned Sequences\n>Q\n%qas\n>T\n%tas\nCoding Sequences\n>Q\n%qcs\n>T\n%tcs\n" >  EXONERATE_OUTPUT.TXT
```

Then, process the output file to extract stop codon and frameshift information:

```
grep -P "Query:|Target:|(Query range:)|(Target range:)|(^[A-Z]\tTAA\t)|(\sframeshift\s)|(^[A-Z]\tTAG\t)|(^[A-Z]\tTGA\t)|(^\*)"  EXONERATE_OUTPUT.TXT   > EXONERATE_ERRORS.TXT
```

Also use this perl tool (https://github.com/EVidenceModeler/EVidenceModeler/blob/master/EvmUtils/misc/Exonerate_to_evm_gff3.pl) by EVM to convert exonerate output file to GFF3 format file.

```
perl Exonerate_to_evm_gff3.pl EXONERATE_OUTPUT.TXT > EXONERATE_OUTPUT.GFF3
```

## Usage

Finally, use all these input files to create a table of possible pseudogenization events:

```
perl tabulate_stops_frameshifts.pl QUERY_CDS_SEQUENCES.fasta EXONERATE_OUTPUT.GFF3 EXONERATE_ERRORS.TXT > STOPS_FRAMESHIFTS.TXT
```
