HistoMapTx Documentation
========================

.. image:: ../figures/logo.png
   :align: right
   :width: 150px

**Spatial Transcriptomics Analysis with Histological Context**

HistoMap is a Python library for analyzing and visualizing histological annotations alongside 
spatially resolved transcriptomics data (Visium). It provides tools for processing, analyzing, 
and visualizing GeoJSON-based tissue annotations with spatial transcriptomics spot data.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorials/tutorials

Quick Start
----------

.. code-block:: python

   import histomap as hm
   import squidpy as sq

   # Load Visium data
   adata = sq.datasets.visium_fluo_moran_test()
   spatial_data = adata.uns['spatial']

   # Load annotations from a GeoJSON file
   histo = hm.HistoMap("annotations.geojson", spatial_data)

   # Plot the annotations
   histo.plot_annotations()

Installation
-----------

.. code-block:: bash

   pip install histomap

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`