.. raw:: html

    <style> .c1 {color:#3282b8; font-weight:bold} </style>

.. role:: c1


overview
=========

**htseq-clip**

htseq-clip is a toolset designed for the processing and analysis of eCLIP/iCLIP dataset.
This package is designed primarily to do the following operations:

:ref:`Annotation <AnnotationOverview>`
***************************************

A suite of functions to process and flatten genome annotation file. 

:c1:`annotation`

:ref:`annotation function <annotation>` takes as input a GFF formatted genome annotation file and converts the annotations from GFF format to bed format.
For an example, this function converts the following GFF annotation

.. _GFFTable:

.. list-table::
   
    * - chr1
      - HAVANA
      - exon
      - 11869
      - 12227
      - .
      - \+
      - .
      - ID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_status=KNOWN;transcript_name=DDX11L1-002;exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1;tag=basic    

and converts this entry into the following BED6 format

.. _BEDTable:

.. list-table::
    :header-rows: 1
    
    * - chromosome
      - start
      - end
      - name
      - score
      - strand
    * - chr1
      - 11868
      - 12227
      - ENSG00000223972.5@DDX11L1@transcribed_unprocessed_pseudogene@exon@1/4@ENSG00000223972.5:exon0001
      - 0
      - \+
Various attributes in the name column in this BED entry is seperated by ``@`` and the
order is given below

.. _AttibTable:

.. list-table::
    :widths: 3,10
    :header-rows: 1
    

    * - atrribute
      - attribute description 
    * - ENSG00000223972.5
      - gene id
    * - DDX11L1
      - gene name
    * - transcribed_unprocessed_pseudogene
      - gene type
    * - exon
      - gene feature (exon, intron, CDS,...)
    * - 1/4
      - 1st exon out of a total of 4 exons of this gene
    * - ENSG00000223972.5:exon0001
      - unique id, merging gene id feature and feature number

:c1:`createSlidingWindows`

:ref:`createSlidingWindows function <createSlidingWindows>` takes as input a flattened annotation BED file
created by the annotation function and splits each individual BED entries into overlapping windows. 
``--windowSize`` parameter controls the size of each window and ``--windowStep`` controls the overlap 
of each neighboring windows from the same feature

Continuing with the example entry above, the first 5 sliding windows generated from the
:ref:`BED6 flattened entry <BEDTable>` are given below:

.. _SWTable:

.. list-table::
    :header-rows: 1
        
    * - chromosome
      - start
      - end
      - name
      - score
      - strand
    * - chr1
      - 11868
      - 11918
      - ENSG00000223972.5@DDX11L1@transcribed_unprocessed_pseudogene@exon@1/4@ENSG00000223972.5:exon0001W00001@1
      - 0
      - \+
    * - chr1
      - 11888
      - 11938
      - ENSG00000223972.5@DDX11L1@transcribed_unprocessed_pseudogene@exon@1/4@ENSG00000223972.5:exon0001W00002@2
      - 0
      - \+
    * - chr1
      - 11908
      - 11958
      - ENSG00000223972.5@DDX11L1@transcribed_unprocessed_pseudogene@exon@1/4@ENSG00000223972.5:exon0001W00003@3
      - 0
      - \+
    * - chr1
      - 11928
      - 11978
      - ENSG00000223972.5@DDX11L1@transcribed_unprocessed_pseudogene@exon@1/4@ENSG00000223972.5:exon0001W00004@4
      - 0
      - \+
    * - chr1
      - 11948
      - 11998
      - ENSG00000223972.5@DDX11L1@transcribed_unprocessed_pseudogene@exon@1/4@ENSG00000223972.5:exon0001W00005@5
      - 0
      - \+

Each sliding window listed here is 50bp long, as default value for ``--windowSize`` argument is ``50``  and the difference between
start positions of each is 20bp, as the default value for ``--windowStep`` argument is ``20`` 

Following the convention in :ref:`flattened annotation <BEDTable>` the attributes in sliding windows name column are also seperated by ``@`` 
and the first 5 attributes in the name column here are exactly the same as that of :ref:`flattened annotation name column <AttibTable>`
An example is given below

.. _SWAttibTable:

.. list-table::
    :header-rows: 1

    * - atrribute
      - attribute description
      - Found in :ref:`flattend name attribute <AttibTable>`
    * - ENSG00000223972.5
      - gene id
      - Yes
    * - DDX11L1
      - gene name
      - Yes
    * - transcribed_unprocessed_pseudogene
      - gene type
      - Yes
    * - exon
      - gene feature (exon, intron, CDS,...)
      - Yes
    * - 1/4
      - 1st exon out of a total of 4 exons of this gene
      - Yes
    * - ENSG00000223972.5:exon0001W00001
      - unique id, merging gene id feature, feature number and window number (W stands for window)
      - No
    * - 1
      - 1st window of this feature 
      - No
 
.. Note:: There will be zero overlap between neighboring windows from two separate gene features

:ref:`Extraction <ExtractionOverview>`
**************************************
Extract and process crosslink sites from alignment file.

:c1:`extract`

:ref:`extract function <extract>` takes as input an alignment file (.bam) and extracts and 
writes either start, insertion, deletion, middle or end site into a BED6 formatted file.
The argument ``--site``  determines crosslink site choice.

.. _AlignTable1:

.. list-table::

  * - HWI-EAS350_213:5:47:1250:7471_4969020:GAAGTCC
    - 0
    - chr1
    - 10063051
    - 255
    - 47M
    - \*
    - 0
    - 0
    - GATGTGTCGGGTACTTGGGCATGAGAGTGAGCAGAGGGAGGAGCTAA
    - \`aa\^a\^aba\]NZY_\^aa]YY\\aV\^P\^PYQ\\RNWLY\`ZMV\[OYNPXT
    - NH:i:1
    - HI:i:1
    - AS:i:46
    - nM:i:0
    - YB:i:3

**start site entry**

.. _StartTable:

.. list-table::

  * - chromosome
    - start
    - end
    - name
    - score
    - strand
  * - chr1
    - 10063050
    - 10063051
    - HWI-EAS350_213:5:47:1250:7471_4969020:GAAGTCC|47
    - 3
    - \+

**middle site entry**

.. _MiddleTable:

.. list-table::

  * - chromosome
    - start
    - end
    - name
    - score
    - strand
  * - chr1
    - 10063074
    - 10063075
    - HWI-EAS350_213:5:47:1250:7471_4969020:GAAGTCC|47
    - 3
    - \+

**end site entry**

.. _EndTable:

.. list-table::

  * - chromosome
    - start
    - end
    - name
    - score
    - strand
  * - chr1
    - 10063097
    - 10063098
    - HWI-EAS350_213:5:47:1250:7471_4969020:GAAGTCC|47
    - 3
    - \+


:ref:`Count <CountOverview>`
****************************
Calculate the number of extracted crosslink sites per given gene annotation feature.

.. figure:: htseq-clip.png
   :width: 75% 