# pegIT

#### About:
pegIT is a python / django tool for quick design of prime editing guide RNAs,
with support for a wide range of targets.

#### Prerequisites:
- [Python](https://www.python.org/download/) - 3.8
- Python libraries are listed in requirements.txt


- [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml "Bowtie")
- [twoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/ "twoBitToFa")
- [primer3](http://primer3.sourceforge.net "primer3")

- For running the web-tool, [celery](https://docs.celeryproject.org), requires a message broker, see [here](https://docs.celeryproject.org/en/stable/getting-started/first-steps-with-celery.html#first-steps) for more info.

make a copy of `CONF.yaml` called `CONF_local.yaml` and update parameters as needed:
- `DESIGN_TWO_BIT_TO_FA_PATH`: path to twoBitToFa
- `DESIGN_BOWTIE_PATH`: path to bowtie
- `DESIGN_PRIMER3_PATH`: path to primer3
- `DESIGN_PRIMER3_CONFIG_PATH`: path to primer3_config
- `DESIGN_ASSEMBLIES_FOLDER`: path to folder containing *.2bit genomes
- `DESIGN_BOWTIE_GENOMES_FOLDER`: path to folder containing bowtie indices. 
- `DESIGN_CODON_USAGE_FOLDER`: path to folder containing codon usage json files
- `DESIGN_OUTPUT_FOLDER`: path to folder to save results to
- `DESIGN_ENTREZ_EMAIL`: email used when searching clinvar database
- `DESIGN_CLINVAR_ORGANISM`: assembly name of the organism that corresponds to clinvar, eg. GCF_000001405.38_GRCh38.p12_genomic

- `DATABASES` can be updated to specifiy alternative database to store data in
- `CELERY_BROKER_URL` should be updated to match with the chosen celery broker.

In order to run, pegIT requires organisms to be added to the database. To add an organism to the database, some files are required:
- A *.2bit genome file
- A GFF3 file containing annotations to be added
- a bowtie index

Before running pegIT for the first time, the database must be instantiated. This can be done by running `python manage.py migrate`

Use [faToTwoBit](http://hgdownload.soe.ucsc.edu/admin/exe/ "faToTwoBit") and [bowtie](http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer) to convert a FASTA file to .2bit and generate bowtie index.

Move the 2bit file to the 2bit folder `DESIGN_ASSEMBLIES_FOLDER`. Move the bowtie indices to the `DESIGN_BOWTIE_GENOMES_FOLDER` should be in a subfolder with same name as the index, eg. `<DESIGN_BOWTIE_GENOMES_FOLDER>homo sapiens/homo sapiens*ebwt`

#### Running pegIT:

pegIT can be run as a command-line tool or a web-tool. To launch the web-tool, run `python manage.py runserver --insecure`, the tool can now be accessed at http://localhost:8000. In web-tool mode, pegit creates oligos in the background, in a separate shell run `python manage.py celery`, to launch celery for processing jobs in the background.


To run the command-line version, run `python pegIT.py`

Available commands in the command-line version
- organisms
    * list - list available organisms
    * add - add organism and organism to the database
    
- design
    * organism - design pegRNAs for desired edits in an organism
    * clinvar - design pegRNAs for introducing or repairing a clinvar variant
    
- utils
    * 2bit - create 2bit file in assembly folder - assumes that faToTwoBit is in same directory as twoBitToFa. Also creates a file with sequence identifiers.
    * bowtie - build bowtie index to bowtie-genomes filder - assumes that bowtie-build is in same directory as bowtie

pegIT requires an edit file containing the desired edits:
- organism mode: a tab-separated file with four columns, without header

        sequence_type	sequence	edit	options
    sequence_type can be on of:
- genomic
- gene
- transcript
- custom

    sequence is the identifier. For genomic the name of a chromosome, for gene the gene_id, for transcript the transcript_id. For custom sequences the sequence field should contain the input sequence and optionally a name in the format \<NAME\>,\<SEQUENCE\>.

    edit is the name of an edit (currently available):
- Insertion
- Deletion
- Substitution
- AminoAcidAlteration
- TagProtein

    options contains a comma separated option string, with<\OPTION\>=<\VALUE/> pairs. eg: insert=A,position=3. specifying repair=true will cause pegIT to design pegRNAs to repair the specified edit.
    
        transcript	HBB-201	AminoAcidAlteration	alteration  =V2A
        gene	HBB	Insertion	insert=ATG,position=200
        genomic	chr3	Deletion	delete=ACAGTCAGTATCAATTCTGGAAGAATTTCCAG,position=46373453

- clinvar mode: a tab-separated file with two columns, without header

        clinvar_id	repair
    
    pegIT can currently only design pegRNAs for clinvar records with a canonical SPDI.
    
    repair can be either true (to design pegRNAs to repair the variant) or false (to design pegRNAs to introduce the variant)
    
    example edit file:
        
        VCV000834692	false
        VCV000834692	true
        
additional parameters (both modes):
- num_pegs - number of pegRNAs to design per edit
- pbs_min_length - minimum length of the pegRNA PBS
- pbs_max_length - maximum legth of the pegRNA PBS
- rt_min_length - minimum length of the homologous region downstream of edit in RT template
- rt_max_length - maximum length of the homologous region downstream of edit in RT template
- primer_min_length - minimum primer length
- primer_opt_length - optimal primer length
- primer_max_length - maximum primer length
- primer_min_tm - minimum primer TM
- primer_opt_tm - optimal primer TM
- primer_max_tm - maximum primer TM
- product_min_size - minimal size of PCR products
- product_max_size - maximal size of PCR products
- spacer_search_range - distance from edit to search for pegRNA spacers
- nicking_range - distance from pegRNA to search for nicking sgRNAs

Example codon usage file:
    
    {
    "TTT": ["F", 0.45], "TTC": ["F", 0.55], "TTA": ["L", 0.07], "TTG": ["L", 0.13],
    "TCT": ["S", 0.18], "TCC": ["S", 0.22], "TCA": ["S", 0.15], "TCG": ["S", 0.06],
    "TAT": ["Y", 0.43], "TAC": ["Y", 0.57], "TAA": ["$", 0.28], "TAG": ["$", 0.2],
    "TGT": ["C", 0.45], "TGC": ["C", 0.55], "TGA": ["$", 0.52], "TGG": ["W", 1.0],
    "CTT": ["L", 0.13], "CTC": ["L", 0.2], "CTA": ["L", 0.07], "CTG": ["L", 0.41],
    "CCT": ["P", 0.28], "CCC": ["P", 0.33], "CCA": ["P", 0.27], "CCG": ["P", 0.11],
    "CAT": ["H", 0.41], "CAC": ["H", 0.59], "CAA": ["Q", 0.25], "CAG": ["Q", 0.75],
    "CGT": ["R", 0.08], "CGC": ["R", 0.19], "CGA": ["R", 0.11], "CGG": ["R", 0.21],
    "ATT": ["I", 0.36], "ATC": ["I", 0.48], "ATA": ["I", 0.16], "ATG": ["M", 1.0],
    "ACT": ["T", 0.24], "ACC": ["T", 0.36], "ACA": ["T", 0.28], "ACG": ["T", 0.12],
    "AAT": ["N", 0.46], "AAC": ["N", 0.54], "AAA": ["K", 0.42], "AAG": ["K", 0.58],
    "AGT": ["S", 0.15], "AGC": ["S", 0.24], "AGA": ["R", 0.2], "AGG": ["R", 0.2],
    "GTT": ["V", 0.18], "GTC": ["V", 0.24], "GTA": ["V", 0.11], "GTG": ["V", 0.47],
    "GCT": ["A", 0.26], "GCC": ["A", 0.4], "GCA": ["A", 0.23], "GCG": ["A", 0.11],
    "GAT": ["D", 0.46], "GAC": ["D", 0.54], "GAA": ["E", 0.42], "GAG": ["E", 0.58],
    "GGT": ["G", 0.16], "GGC": ["G", 0.34], "GGA": ["G", 0.25], "GGG": ["G", 0.25]
    }
