#!/bin/bash

input="toDL.csv"

while IFS= read -r line 
do
	echo $line
	wget --no-check-certificate $line 
 done < $input 


