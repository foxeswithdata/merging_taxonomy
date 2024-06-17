#!/bin/bash

input="dirs_to_create.csv"

while IFS= read -r line 
do
	echo $line
	mkdir -p  $line 
 done < $input 


