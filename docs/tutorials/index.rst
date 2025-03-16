Tutorials
=========

Getting Started with HistoMapTx
-------------------------------

These tutorials guide you through the process of analyzing spatial transcriptomics data with histological annotations using HistoMapTx. They are designed to be followed in sequence, building your skills from basic annotation handling to advanced spatial analysis techniques.

Each tutorial includes practical examples with sample data and explanations of key concepts in spatial transcriptomics analysis.

Basic Concepts
-------------

Before diving into the tutorials, it's helpful to understand these key concepts:

- **Annotations**: Histological regions marked on tissue images (e.g., tumor, stroma, necrosis)
- **Spots**: Spatial locations on a Visium slide where gene expression is measured
- **Annotation Map**: Assignment of spots to specific tissue annotations
- **Overlap**: Quantification of how much each spot intersects with annotations

Tutorial Series
--------------

.. toctree::
   :maxdepth: 1
   :caption: Tutorial Notebooks:

   ../notebooks/1_Basic_annotations.ipynb
   ../notebooks/2_Mapping_annotations.ipynb
   ../notebooks/3_Working_with_annotations.ipynb
   ../notebooks/4_Adding_segmentation.ipynb

Tutorial Descriptions
--------------------

**1. Basic Annotations**
   Learn how to load and visualize histological annotations from sources like QuPath. This tutorial covers tissue detection, visualizing annotations on histology images, and managing annotation properties.
   
**2. Mapping Annotations to Spots**
   Discover how to assign tissue annotations to Visium spots. Topics include calculating overlap between spots and annotations, visualizing annotation coverage, and generating annotation maps for downstream analysis.
   
**3. Working with Annotations**
   Explore advanced spatial analysis techniques using annotation maps. Learn to calculate distances from tissue regions, analyze gene expression gradients relative to boundaries, and identify genes with specific spatial patterns.
   
**4. Adding Cell Segmentation**
   Integrate cell segmentation data with spatial transcriptomics. This tutorial shows how to add cell information to your analysis, visualize cell density across tissue regions, and analyze cellular features in different annotations.

Getting Help
-----------

If you encounter issues while following these tutorials, please:

- Check the :ref:`API documentation <api>` for detailed function descriptions
- Visit our `GitHub repository <https://github.com/Dantferno/HistoMapTx>`_ for the latest updates
- Report bugs or ask questions in the GitHub issues section