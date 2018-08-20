5. Examples
==================


---------------------------------
5.1 TFBS_footprinter Use Examples
---------------------------------
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Running the sample analyses
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Run the sample analysis using a .csv of arguments:
	
	``$ tfbs_footprinter -t PATH_TO/sample_analysis/sample_analysis_list.csv``

- Run the sample analysis using a .txt of Ensembl transcript ids, and minimal arguments:
	
	``$ tfbs_footprinter -t PATH_TO/sample_analysis/sample_ensembl_ids.txt``

- Run the sample analysis using a .txt of Ensembl transcript ids, and all arguments:

	``$ tfbs_footprinter -t PATH_TO/sample_analysis/sample_ensembl_ids.txt -tfs PATH_TO/sample_analysis/sample_jaspar_tf_ids.txt -s homo_sapiens -g mammals -e low -pb 900 -pa 100 -tx 10 -update``


---------------------------------------------------
5.2 TFBS_footprinter Use Examples **Within Docker**
---------------------------------------------------
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Running the sample analyses
^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Within Docker we first need to mount a volume so that the results of the analyis can be viewed on our host computer.  It is recommended that you create an empty directory on your host computer:

	``$ docker run -v /ABSOLUTE_PATH_TO/EMPTY_DIR_ON_HOST:/home/sample_analysis/tfbs_results -it tfbs_footprinting bash``

2. Then we move into the pre-existing sample analysis directory in the Docker container to perform the analysis there so that the results generated there will automatically appear in the designated location on our host computer:

	``$ cd ./sample_analysis``

3. Then we can run the sample analysis in Docker in the same way that we would normally use tfbs_footprinter (above), e.g. using a .csv of arguments:

	``$ tfbs_footprinter -t ./sample_analysis_list.csv``

- Or (again, as above) using a .txt of Ensembl transcript ids, and minimal arguments:

	``$tfbs_footprinter -t ./sample_ensembl_ids.txt``

- Or (again, as above) using a .txt of Ensembl transcript ids, and multiple arguments:

	``$ tfbs_footprinter -t ./sample_ensembl_ids.txt -tfs ./sample_jaspar_tf_ids.txt -s homo_sapiens -g mammals -e low -pb 900 -pa 100 -tx 10 -o ./tfbs_results -update``


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Example using user-defined files/arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Within Docker we first need to mount a volume so that we can load your analysis files from your host computer into docker AND save the results of the analysis on our host computer:

	``$ docker run -v /ABSOLUTE_PATH_TO/DIR_ON_HOST/CONTAINING_ANALYSIS_FILES:/home/analysis_dir -it tfbs_footprinting bash``
2. Then we move into your analysis directory in the Docker container to perform the analysis there so that the results generated there will automatically appear in the designated location on our host computer:

	``$ cd ./analysis_dir``
3. Then we can run the sample analysis in Docker in the same way that we would normally use tfbs_footprinter (above), e.g. using a .csv of arguments:

	``$ tfbs_footprinter -t ./USER_TABLE_OF_ENSEMBL_IDS_AND_ARGS.csv``

- Or (again, as above) using a .txt of Ensembl transcript ids, and minimal arguments:

	``$ tfbs_footprinter -t ./USER_LIST_OF_ENSEMBL_IDS.txt``

- Or (again, as above) using a .txt of Ensembl transcript ids, and multiple arguments:

	``$ tfbs_footprinter -t ./USER_LIST_OF_ENSEMBL_IDS.txt -tfs ./USER_LIST_OF_TF_NAMES.txt -s homo_sapiens -g mammals -e low -pb 900 -pa 100 -tx 10 -o PATH_TO/Results/ -update``


----------------------------------
5.3 Update experimental data files
----------------------------------
This does not need to be run once every few months, check the documentation on `ReadTheDocs <https://tfbs-footprinting.readthedocs.io/en/latest/>`_ for new releases.

``$ tfbs_footprinter -update``