9. Troubleshooting
==================

---------------------
9.1 Log file
---------------------
The first step in troubleshooting any possible errors that might arise is to check the log file.  The log file is created in the directory where you initiated the analysis, and it is named 'TFBS_footprinter.log'.  Many relevant events are logged there instead of output to terminal:

- start time of analysis
- arguments/settings used in analysis
- full queries made to the Ensembl REST system
- warnings from Ensembl REST system regarding exceeding of rate limit
- if given transcript ids are not in the Ensembl system (possibly misformed or deprecated)
- if experimental data used in TFBS prediction has been downloaded
- if there was an error in downloading experimental data
- if results already exist for the current analysis, which are then loaded
- total time of analysis
- if there was an error in retrieving an alignment for the defined region
- if the transcript information file already exists
- if the transcript regulatory information file already exists

---------------------
9.1 Docker
---------------------
The Docker installation will have a default RAM allocation that is too low (~2GB) to run TFBS_footprinting.  This setting should be changed to >6144MB.

In Windows this can be adjusted by navigating through: 

**Docker system tray icon>Settings>Advanced>Memory**.  

After changing this value Docker will restart, which could take several minutes.  If the allocated memory is too low, then the Docker image will terminate and you will see the ``$ killed`` message in the console.




