#!/bin/bash

RUNID=$1

zip -r ${RUNID}.zip ./${RUNID}/model_summaries/ ./${RUNID}/model_predictions/ ./${RUNID}/R_output/ ./${RUNID}/run_gams.R ./${RUNID}/run_gams.R.cmd ./${RUNID}/*.Rdata