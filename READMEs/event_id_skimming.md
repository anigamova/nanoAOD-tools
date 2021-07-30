# Event ID Skimming

If one wishes to derive scaling equations for analysis categories, one needs a reweighted dataset of events which are binned according to the categories. Of course, this could be done by writing an analysis and running over the NanoAOD, binning the events as you go. However, this could prove to be a lengthy task. Thankfully, there exists a better approach!

 An event in a centrally-produced dataset will have a unique ID that is carried across data tiers. Suppose that an analysis has already been performed and an output of that analysis is a list of event IDs and which analysis category that event belongs to, then we have everything we need. We can find the NanoAOD version of whatever dataset was input to the analysis and simply pick out the events we want (based on their ID) and add the category tags as we go along. Then we can reweight this skimmed dataset as usual (see [walkthrough tutorial](walkthrough.md)).

 This task is performed by the [eventID_skim.py](../eventIDSkimming/eventID_skim.py) script. The event IDs are input to the script via a csv file. This file needs at least two columns, named `event` and `category` which correspond to the event IDs and category tags respectively. If wanted, further columns can be included which can be added to the NanoAOD.
 
 The script is run in the following way:
 ```
python eventIDSkimming/eventID_skim.py theDataset eventIDs.csv -o output.root
 ```
 Options include:
 - `--extraBranches` a separated list of branches/columns from the csv file that should be added to the output tree. Using 'auto' adds all branches from the csv file.
 - `--extraCollections` a comma separated list of variables/collections from the nanoAOD to include. All collections needed for reweighting are included by default, this option is for additional collections. 
 - `--test` sets the scrip to skim over only one file and up to 10,000 events

Note that the full name of the dataset need not be used. The [datasetKeys.json](../eventIDSkimming/datasetKeys.json) file exists so shorter names can be used instead. The script will work if you use the shortened or long version of the name. If the input dataset is not found in [datasetKeys.json](../eventIDSkimming/datasetKeys.json) it will assume you have an input a long version and will proceed as normal.