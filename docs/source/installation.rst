requirements and installation
=============================

requirements
************

* Python >= 3.5
* `HTSeq`_  package

.. Note:: `HTSeq`_ uses `pysam`_ package for processing alignment files.  Please consult `HTSeq manual`_ and  `pysam manual`_ for requirements of both packages.

.. _`HTSeq`: https://pypi.org/project/HTSeq/
.. _`pysam`: https://pypi.org/project/pysam/
.. _`HTSeq manual`: https://htseq.readthedocs.io/en/release_0.11.1/install.html
.. _`pysam manual`: https://pysam.readthedocs.io/en/latest/installation.html

installation
************

quick installation
------------------

If a user has a local python environment with all the dependencies for HTSeq and pysam installed,
then htseq-clip can be installed as:

.. code-block:: sh 

    $ pip install htseq-clip

conda environment
-----------------

We strongly encourage the use of `conda`_ package management system for multiple Python versions/various incompatible package installations.
Please install conda on your computer following `the guidelines`_. Once conda installation is successful, create a new ``htseq-clip`` enviroment as:

.. _`conda`: https://docs.conda.io/en/latest/
.. _`the guidelines`: https://docs.anaconda.com/anaconda/install/

.. code-block:: sh 

    (base) $ conda create -n htseq-clip

and activate the environment:

.. code-block:: sh 

    (base) $ conda activate htseq-clip

now install the dependencies:

.. code-block:: sh 

    (htseq-clip) $ conda install -c bioconda pysam
    ....
    (htseq-clip) $ conda install -c bioconda htseq

now htseq-clip can be installed in this environment as:

.. code-block:: sh 

    (htseq-clip) $ pip install htseq-clip
