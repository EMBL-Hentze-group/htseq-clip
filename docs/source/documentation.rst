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

Annotation
**********

  * :c1:`annotation`

    Flattens a given annotation file in GFF format to BED6 format 
    
    **Arguments**

    * ``-g/--gff``       GFF formatted annotation file, supports .gz files
    * ``-u/--geneid``    Gene id attribute in GFF file 
    * ``-n/--genename``  Gene name attribute in GFF file
    * ``-t/--genetype``  Gene type attribute in GFF file
    * ``--splitExons``   This flag splits exons into components such as 5' UTR, CDS and 3' UTR
    * ``--unsorted``     Use this flag if the GFF file is unsorted
    * ``-o/--output``    Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console   

    The default values for ``--geneid``, ``--genename`` and ``--genetype`` arguments follow 
    `gencode GFF format`_
    
    .. _`gencode GFF format`: https://www.gencodegenes.org/pages/data_format.html

    **Usage**
    
    .. code-block:: sh    
      
      $ htseq-clip annotation -h  

  * :c1:`createSlidingWindows`

    Create sliding windows from the flattened annotation file

    **Arguments**

    * ``-i/--input``  Flattened annoation file, see 
    * ``-w/--windowSize``  Window size in number of base pairs for the sliding window
    * ``-s/--windowStep``  Window step size for sliding window
    * ``-o/--output``    Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console


    **Usage**

    .. code-block:: sh    
      
      $ htseq-clip createSlidingWindows -h
    
  * :c1:`mapToId`

    Extract "name" column from the annotation file and map the entries to unique id 
    and print out in tab separated format

    **Arguments**

    * ``-a/--annotation``  Flattened annotation file from --- or sliding window file from ---
    * ``-o/--output``    Output file name. If the file name is given with .gz suffix, it is gzipped. If no file name is given, output is print to console

    
    **Usage**
    
    .. code-block:: sh    
      
      $ htseq-clip mapToId -h

Extraction
**********

  * :c1:`extract`

    Extract crosslink sites, insertions or deletions

    .. code-block:: sh    
      
      $ htseq-clip extract -h

Counting
********
  
  * :c1:`count`

    Counts the number of crosslink/deletion/insertion sites

    .. code-block:: sh    
      
      $ htseq-clip count -h

Helpers
*******
  
  * :c1:`createMatrix`
    
    Create R friendly output matrix file from count function output files

    .. code-block:: sh    
      
      $ htseq-clip createMatrix -h