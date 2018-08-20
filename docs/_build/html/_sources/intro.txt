.. _intro:

1. Introduction
==================
----------
1.1 Basics
----------
**Primary Usage:** Identification of cis-regulatory elements initially identified by matrix scoring and then additionally scored on 7 other relevant contextual datapoints.  Based on analysis of protein-coding transcripts in the Ensembl database.

.. image:: https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/tfbs_logo.png
	:alt: logo

- This work is a derivative of `"Transcription factors" <https://commons.wikimedia.org/wiki/File:Transcription_Factors.svg>`_ by `kelvin13 <https://commons.wikimedia.org/wiki/User:Kelvin13>`_, used under `CC BY 3.0 <https://creativecommons.org/licenses/by/3.0/>`_.


.. note::
	TFBS_footprinting is now available in a `Docker image <https://hub.docker.com/r/thirtysix/tfbs_footprinting/>`_.


---------------------
1.2 Experimental data
---------------------
The TFBS footprinting method computationally predicts transcription factor binding sites (TFBSs) in a target species (e.g. homo sapiens) using 575 position weight matrices (PWMs) based on binding data from the JASPAR database.  Additional experimental data from a variety of sources is used to support or detract from these predictions:

* DNA sequence conservation in homologous mammal species sequences
* proximity to CAGE-supported transcription start sites (TSSs)
* correlation of expression between target gene and predicted transcription factor (TF) across 1800+ samples
* proximity to ChIP-Seq determined TFBSs (GTRD project)
* proximity to qualitative trait loci (eQTLs) affecting expression of the target gene (GTEX project)
* proximity to CpGs
* proximity to ATAC-Seq peaks (ENCODE project)



--------------------------
1.3 Ensembl transcript-ids
--------------------------
In order to incorporate several of the experimental datapoints, and for ease of use, the Ensembl transcript ID was chosen as the primary identifier around which TFBS_footprinting analyses would revolve.  This allows users to predict TFBSs in the promoters any of 1-80,000 human protein coding transcripts in the Ensembl database, and incorporate expression data from the FANTOM project, and eQTLS from the GTEx project, which has been assigned to these IDs.  TFBS predictions can also be made for 87 unique non-human species (including model organisms such as mouse and zebrafish), present in the following groups:

- 70 Eutherian mammals
- 24 Primates
- 11 Fish
- 7 Sauropsids

