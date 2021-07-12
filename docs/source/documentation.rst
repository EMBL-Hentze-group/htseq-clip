.. raw:: html

    <style> .c1 {color:#3282b8; font-weight:bold} </style>

.. role:: c1

documentation
=============

After successful installation of the package use 

.. code-block:: sh

   $ htseq-clip -h

for a brief description of the functions available in htseq-clip. 
The available functions can be categorized into 4 different classes given below.


.. _AnnotationOverview:

Annotation
**********

.. _annotation:

:c1:`annotation`
-----------------

  Flattens a given annotation file in GFF format to BED6 format 
  
  **Arguments**

  * ``-g/--gff``       GFF formatted annotation file, supports .gz files
  * ``-u/--geneid``    Gene id attribute in GFF file (default: gene_id)
  * ``-n/--genename``  Gene name attribute in GFF file (default: gene_name)
  * ``-t/--genetype``  Gene type attribute in GFF file (default: gene_type)
  * ``--splitExons``   This flag splits exons into components such as 5' UTR, CDS and 3' UTR
  * ``--unsorted``     Use this flag if the GFF file is unsorted
  * ``-o/--output``    Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console   

  .. Note:: The default values for ``--geneid``, ``--genename`` and ``--genetype`` arguments follow  `gencode GFF format`_

  .. _`gencode GFF format`: https://www.gencodegenes.org/pages/data_format.html


  **Usage**
  
  .. code-block:: sh    
    
    $ htseq-clip annotation -h  

.. _createSlidingWindows:

:c1:`createSlidingWindows`
---------------------------

  Create sliding windows from the flattened annotation file

  **Arguments**

  * ``-i/--input``  Flattened annoation file, see annotation_
  * ``-w/--windowSize``  Window size in number of base pairs for the sliding window (default: 50)
  * ``-s/--windowStep``  Window step size for sliding window (default: 20)
  * ``-o/--output``    Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console


  **Usage**

  .. code-block:: sh    
    
    $ htseq-clip createSlidingWindows -h
  
.. _mapToId:

:c1:`mapToId`
-------------

  Extract "name" column from the annotation file and map the entries to unique id 
  and print out in tab separated format

  **Arguments**

  * ``-a/--annotation``  Flattened annotation file from annotation_ or sliding window file from createSlidingWindows_
  * ``-o/--output``    Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console


  **Usage**

  .. code-block:: sh    
    
    $ htseq-clip mapToId -h

.. _ExtractionOverview:

Extraction
**********

.. _extract:

:c1:`extract`
-------------

  Extract crosslink sites, insertions or deletions

  **Arguments**

  * ``-i/--input`` Input .bam file. Input bam file must be co-ordinate sorted and indexed
  * ``-e/--mate`` for paired end sequencing, select the read/mate to extract the crosslink sites from, accepted choices: ``1, 2``

    * ``1`` use the first mate in pair
    * ``2`` use the second mate in pair
  * ``-s/--site`` Crosslink site choices, accepted choices: ``s, i, d, m, e`` (default: e)
    
    * ``s`` startsite, 
    * ``i`` insertion site 
    * ``d`` deletion site 
    * ``m`` middle site 
    * ``e`` end site 
  
  * ``-g/--offset`` Number of nucleotides to offset for crosslink sites (default: 0)
  * ``--ignore`` Use this flag to ignore crosslink sites outside of genome annotations
  * ``-q/--minAlignmentQuality`` Minimum alignment quality (default: 10)
  * ``-m/--minReadLength`` Minimum read length (default: 0)
  * ``-x/--maxReadLength`` Maximum read length (default: 500)
  * ``-l/--maxReadInterval`` Maximum read interval length (default: 10000)
  * ``--primary`` Use this flag consider only primary alignments of multimapped reads
  * ``-c/--cores`` Number of cores to use for alignment parsing (default: 5)
  * ``-t/--tmp`` Path to create and store temp files (default behavior: use parent folder from "--output" parameter)
  * ``-o/--output`` Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console
  
  **Usage**

  .. code-block:: sh    
    
    $ htseq-clip extract -h
  
  
  .. Note:: To extract ``1``st offset position of second mate (``2``) start site (``s``) in eCLIP, use: ``--mate 2 --site s --offset -1``


.. _CountOverview:

Counting
********
  
.. _count:

:c1:`count`
-----------

  Counts the number of crosslink/deletion/insertion sites

  **Arguments**

  * ``-i/--input`` Extracted crosslink sites, see extract_
  * ``-a/--ann`` Flattened annotation file, see annotation_ OR sliding windows file, see createSlidingWindows_
  * ``--unstranded`` crosslink site counting is strand specific by default. Use this flag for non strand specific crosslink site counting
  * ``-o/--output`` Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console

  **Usage**

  .. code-block:: sh    
    
    $ htseq-clip count -h

Helpers
*******
  
.. _createMatrix:

:c1:`createMatrix`
-------------------
    
  Create `R`_ friendly output matrix file from count function output files

  .. _`R`: https://www.r-project.org/

  **Arguments**

  * ``-i/--inputFolder`` Folder name with output files from count function, see count_
  * ``-b/--prefix`` Use files only with this given file name prefix (default: None)
  * ``-e/--postfix`` Use files only with this given file name postfix (default: None)
  * ``-o/--output`` Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console

  .. Warning:: either ``--prefix`` or ``--postfix`` argument must be given

  **Usage**

  .. code-block:: sh    
    
    $ htseq-clip createMatrix -h