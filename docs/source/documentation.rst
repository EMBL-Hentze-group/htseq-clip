.. raw:: html

    <style> .c1 {color:#3282b8; font-weight:bold} </style>

.. role:: c1

documentation
=============

After successful installation of the package use 

.. code-block:: sh

   $ htseq-clip -h

for a brief description of the functions available in htseq-clip. 
The available functions can be categorized into 4 different classes:

**Annotation**
**************

  * :c1:`annotation`

    Flattens a given annotation file in GFF format to BED6 format 

    .. code-block:: sh    
      
      $ htseq-clip annotation -h  

  * :c1:`createSlidingWindows`

    Create sliding windows from the flattened annotation file   

    .. code-block:: sh    
      
      $ htseq-clip createSlidingWindows -h
    
  * :c1:`mapToId`

    Extract "name" column from the annotation file and map the entries to unique id 
    and print out in tab separated format

    .. code-block:: sh    
      
      $ htseq-clip mapToId -h

**Extraction**
**************

  * :c1:`extract`

    Extract crosslink sites, insertions or deletions

    .. code-block:: sh    
      
      $ htseq-clip extract -h

**Counting**
************
  
  * :c1:`count`

    Counts the number of crosslink/deletion/insertion sites

    .. code-block:: sh    
      
      $ htseq-clip count -h

**Helpers**
***********
  
  * :c1:`createMatrix`
    
    Create R friendly output matrix file from count function output files

    .. code-block:: sh    
      
      $ htseq-clip createMatrix -h