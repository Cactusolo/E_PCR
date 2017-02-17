
           Electronic PCR commandline tools: operating instructions

                                                          Version: 2.3.12
     _________________________________________________________________

   Use e-PCR to map sequences using STS database

   Use re-PCR to map STSes or short primers in sequence database

   Use famap and fahash to prepare sequence database for re-PCR searches.
     _________________________________________________________________

Forward e-PCR

Example

work> e-PCR -w9 -f 1 -m100 mystsdb.sts D=100-400 myfastafile.fa N=1 G=1 T=3

Synopsis


e-PCR [-hV] [posix-options] stsfile [fasta ...] [compat-options]
where posix-options are:
        -m ##   Margin (default 50)
        -w ##   Wordsize  (default 7)
        -n ##   Max mismatches allowed (default 0)
        -g ##   Max indels allowed (default 0)
        -f ##   Use ## discontiguos words
        -o ##   Set output file
        -t ##   Set output format:
                1 - classic, range (pos1..pos2)
                2 - classic, midpoint
                3 - tabular
                4 - tabular with alignment in comments (slow)
        -d ##-## Set default sts size
        -p +-   Turn hits postprocess on/off
        -v +-   Verbose on/Off
        -a a|f  Use presize alignmens (only if gaps>0), slow
                 a - Allways or f - as Fallback
        -x +-   Use 5'-end lowercase masking of primers (default -)
        -u +-   Uppercase all primers (default -)
and compat-options (duplicate posix-options) are:
        M=##    Margin (default 50)
        W=##    Wordsize  (default 7)
        N=##    Max mismatches allowed (default 0)
        G=##    Max indels allowed (default 0)
        F=##    Use ## discontinuos words
        O=##    Set output file to ##
        T=##    Set output format (1..4)
        D=##-## Set default sts size
        P=+-    Postprocess hits on/off
        V=+-    Verbose on/Off
        A=a|f   Use presize alignmens (only if gaps>0), slow
                 a - Allways or f - as Fallback
        X=+-    Use 5'-end lowercase masking of primers (default -)
        U=+-    Uppercase all primers (default -)
        -mid    Same as T=2

Description

   e-PCR parses stsfile in unists format, then reads nucleotide sequence
   data in FASTA format from files listed in commandline if any, or from
   stdin otherwise. For input sequences e-PCR finds matches and prints
   output in one of three formats.

Options

   Two sets of options are used: POSIX-compatible and old-style provided
   for compatibility with previous versions of e-PCR.

   Posix-style options can appear only before first parameter not
   starting with '-'. Argument '--' explicitely stops parsing arguments
   as posix options.

   Compatibility options can appear anywhere in commandline. '-mid' can
   appear anywhere and do not stop posix options recognision.

General options

   -V
          Print version, exit after parsing commandline

   -h
          Print help, exit after parsing commandline

Hash building options

   -w wordsize | W=wordsize
          Set word size for primers hash (nucleotide positions). Longer
          word size decreases hash collision rate, but increases memory
          usage. Also no mismatches are allowed within word size near
          "inner" boundary of primers unless one uses discontiguous
          words, and no gaps are ever allowed in that region.

   -f wordcnt | W=wordcnt
          Set discontiguous word count for primers hash (1 means "use
          contiguous words"). Discontiguous words increase number of hash
          tables and decrease "effective" word size (thus increasing hash
          collision rate), so make search significantly slower, but
          increase sencitivity by allowing mismatches within word size.
          Reasonable values are 1 (contiguous words) and 3.

   -d lo-hi | D=lo-hi
          Set ddefault STS size range - values used for STSs that have no
          size associated in file.

Hit quality options

   -m margin | M=margin
          Set maximal allowed deviation of hit product size from expected
          STS size.

   -n mism | N=mism
          Set maximal number of mismatches allowed in primer-to-sequence
          alignment (per primer!).

   -g mism | G=mism
          Set maximal number of gaps allowed in primer-to-sequence
          alignment (per primer!).

Alignment algorithms options

   -a a|f | A=a|f
          Use NW algorithm to align primers to sequence: a - always, f -
          as fallback if fast algorithm gives no hit at this position.

   -x +|- | X=+|-
          Turn on/off recognising of lowercase characters at 5'-ends of
          primers as nucleotides that don't need to be aligned to
          sequence (floppy tails).

   -u +|- | U=+|-
          Uppercase primers. To use with files prepared for ``-x=+''
          mode, but requiring full primer alignment.

   If STS file contains primers with lowercase charactars, you have to
   use either -x+ or -u+ flag.

Report options

   -o output | O=output
          Set output file.

   -t 1|2|3|4 | T=1|2|3|4
          Set output format.

   -p +|- | P=+|-
          Set hit grouping on/off: when using discontiguous words and
          gaps, some hits may be reported multiple times with little
          different quality. This option controls reporting only best hit
          of group of overlapping hits. Default depends on F and G
          values.

   -v +|- | V=+|-
          Report sequence ids to stderr on/off.

Ouput formats

   1: Traditional: reports whitespace-separated

          + Sequence FASTA identifier
          + POS1..POS2 -- start and end positions of hit (includes length
            floppy tail)
          + STS identifier (col. 1 from STS file)
          + STS description (columns 5..last from STS file)

          In this format product size equals to POS2-POS1+1

   2: Traditional midpoint: reports whitespace-separated

          + Sequence FASTA identifier
          + POS -- middle point position of hit
          + STS identifier (col. 1 from STS file)
          + STS description (columns 5..last from STS file)

   3: Tab-separated detailed

          + Sequence FASTA identifier
          + STS identifier (col. 1 from STS file)
          + +|- -- strand of hit (order of primers in hit)
          + POS1 -- start position of hit (does not include floppy tail
            if any)
          + POS2 -- end position of hit (does not include floppy tail)
          + SIZE/MIN..MAX -- observed size of hit/expected size range of
            STS
          + MISM -- Total number of mismatches for two primers
          + GAPS -- Total number of gaps for two primers
          + STS description (columns 5..last from STS file)

          In this format product size may be greater then POS2-POS1+1 for
          probes with floppy tails

   4: Tab-separated detailed with alignment
          Is same as format 3, but also containing visualisations of
          alignments in comment lines (lines starting with ``#'')

Exit codes

   Zero on success, nonzero on fail
     _________________________________________________________________

Reverse e-PCR

Example

work> famap -tN -b genome.famap org/chr_*.fa
work> fahash -b genome.hash -w 12 -f3 ${PWD}/genome.famap
work> re-PCR -s genome.hash -n1 -g1 ACTATTGATGATGA AGGTAGATGTTTTT 120-200

Synopsis


famap [-hV]
famap -b mmapped-file [-t cvt] [fasta-file ...]
famap -d mmapped-file [ord ...]
famap -l mmapped-file [ord ...]
where cvt is one of: off n N nx NX

fahash [-hV]
fahash -b hash-file [build-options] mmapped-file ...
fahash -T hash-file [-o output]

where:
        -b hash-file    Build hash tables (hash-file) from sequence files,
        -T hash-file    Print word usage statistics for hash-file
        -o outfile      Set output file name for -T

build-options:
        -w wordsize     Set word size when building hash tables
        -f period       Set discontiguity when building hash tables
        -k              Skip repeats when building indexfile
        -F min,max      Set watermarks for fragment size (in Mb) for -v1
        -v 1|2          Build file of format version 1 or 2
        -c cachesize    Use cache size cachesize (for -v2)

re-PCR [-hV]
re-PCR -p hash-file [-g gaps] [-n mism] [primer ...]
re-PCR -P hash-file [-g gaps] [-n mism] [primer-file ...]
re-PCR -s hash-file [search-options] [-O output] [left right lo hi [...]]
re-PCR -S hash-file [search-options] [-O output] [-C bcnt] [stsfile ...]

where:
        -p hash-file    Perform primer lookup using hash-file
        -P hash-file    Perform primer lookup using hash-file
        -s hash-file    Perform STS lookup using hash-file, STSs in cmdline
        -S hash-file    Perform STS lookup using hash-file, STSs in file


search-options:
        -n mism         Set max allowed mismatches per primer for lookup
        -g gaps         Set max allowed indels per primer for lookup
        -m margin       Set variability for STS size for lookup
        -d min-max      Set default STS size (for STSs without size set)
        -r +|-          Enable/disable reverse STS lookup
        -O +|-          Enable/disable syscall optimisation

        -C batchcnt     Set number of STSes per batch
        -o outfile      Set output file name

Description

   Reverse e-PCR (re-PCR) performs STS or primer lookup against sequence
   database. Two files are required for database: mmapped-file with
   sequence data in fast random-accessible format and hash-file, that
   keeps precalculated positions of all words of sequence database

   Use famap to build mmapped-file from FASTA files.

   Use fahash to build hash-file, and output word usage statistics.

   Use re-PCR to perform STS and primer searches.

   Discontiguous words are supported by re-PCR as well as contiguous.

Options

Common options

   -V
          Print version, exit after parsing commandline

   -h
          Print help, exit after parsing commandline

famap options

   -b mmapped-file
          Build famap-file from input fasta file(s). If no fasta files
          are set in commandline, use stdin as input.

   -d mmapped-file
          Dump famap-file contents in fasta format. If ord number(s) are
          set, print only sequences with given ordinals.

   -l mmapped-file
          List fama-file sequence identifiers. If ord number(s) are set,
          print only sequences with given ordinals.

   -t cvt-table
          Use compiled-in table to convert input.

        n
                Nucleotides. Allowed characters are [actgACTGnN]. Other
                letters are converted to n or N. Rest of symbols are
                ignored. Case is preserved.

        nx
                Nucleotides with extended ambiquity codes iupac_na,
                lowercase are allowed. Other letters are converted to n
                or N. Rest of symbols are ignored. Case is preserved.

        N
                Nucleotides. Allowed characters are [ACTGN]. [actgn] are
                converted to uppercase. Other letters are converted to N.
                Rest of symbols are ignored.

        NX
                Nucleotides with extended ambiquity codes iupac_na,
                lowercase are converted to uppercase. Other letters are
                converted to N. Rest of symbols are ignored.

Fahash

   -b hash-file
          Build hash-file for mmapped-file(s).

   -T hash-file
          Dump word usage statustics for hash-file.

   -v version
          Build hash-file of version 1 or 2 (2 is default).

   -w wordsize
          Build hash-file for word wordsize nucleotides long.

   -f wordcnt
          Build hash-file for wordcnt discontiguous words. 1 stands for
          contiguous words.

   -F min,max
          Use memory watermarks (Mbytes) for hash table size (for -v 1).

   -c cachesize
          Set cache size for -v 2.

   -o output-file
          Use output-file for output result of -T.

Commands

   -p hash-file
          Perform lookup for primers given in commandline.

   -s hash-file
          Perform lookup for STSes given in commandline.

   -S hash-file
          Perform lookup for STSes taken from unists file(s) given in
          commandline.

Search options

   -n mism
          Number of mismatches allowed per primer.

   -g gaps
          Number of gaps allowed per primer.

   -m margin
          Maximal deviation of observed product size to expected STS
          size.

   -d lo-hi
          Set ddefault STS size range - values used for STSs that have no
          size associated in file.

   -r +|-
          Enable|disable flipped STS lookup (default is "enabled").

   -O +|-
          Enable|disable syscall optimisation. Since lookup is i/o
          expensive, enabling this parameter may improve search
          performance diskwise. On the other hand, it takes significantly
          more memory and CPU.

   -C batchcount
          How many STSs from input file to look at one pass. May effect
          on performance, especialy when used with -O +.

   -o output-file
          Use output-file for output.

Output format

   Is tab-separated file with following fields:

For primer lookup

     * Primer ID
     * Sequence ID
     * Strand
     * Hit start
     * Hit end
     * Mismatches
     * Gaps
     * Size

For STS lookup

     * STS ID
     * Sequence ID
     * Strand
     * Hit start
     * Hit end
     * Mismatches
     * Gaps
     * Observed Size/Expected size range

Exit codes

   Zero on success, non-zero on errors

Bugs and features

     * Mmapped-file path is hardcoded to hash-file as it is in
       commandline when hash-file is being built, which means that when
       one performs searches mmapped-file should be accessible with same
       name from current directory, as it is hardcoded.
     * Mmapped-file is a proprietary format, that could be substituted
       with megablast database format, but is not (yet?) for performance
       reasons.
     * If sequence sizes are large, it may be tricky to create database
       with discontiguous words because of memory usage requirements.
       Changing parameter -F (for -v 1) or -c (for -v 2) may help.
     _________________________________________________________________

File formats

   STS database
          Is single-tab (i.e. two tabs in a row mean "empty field")
          separated file with following fields:

          + STS id (required).
          + First (left) primer (required).
          + Second (right) primer (required).
          + Product size (optional): can be number for strict size, or
            two numbers separated by dash for size range.
          + Additional info, that can be used by applications (optional).

          Primers should be in iupac_na encoding, everything that is not
          ACTG or actg is translated to N or n. Primers sequences should
          be uppercase, unless you want to use file with e-PCR -x+ flag -
          then several first nucleotides of primers may be
          lowercase-masked. If primers are not fully uppercase and you
          don't use -x+ flag, you have to use -u+ flag with e-PCR.

   Primers file
          Is single-tab (i.e. two tabs in a row mean "empty field")
          separated file with following fields:

          + Primer id (required).
          + Primer sequence.
     _________________________________________________________________
