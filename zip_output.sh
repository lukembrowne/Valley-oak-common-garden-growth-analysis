#!/bin/bash

RUNID=$1

zip -r ${RUNID}.zip ./${RUNID}/model_summaries/ ./${RUNID}/R_output/ ./${RUNID}/run_gams.R ./${RUNID}/run_gams.R.cmd ./${RUNID}/*.Rdata