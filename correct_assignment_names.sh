#!/bin/bash

f() { array=("${BASH_ARGV[@]}"); }

for student_dir in *; do

	if [[ -d "$student_dir" ]]; then

		# replace ' ' to '_' in directory names
		underscored_student_dir=`echo ${student_dir:0:${#student_dir}-1} | tr ' ' '_'`

		IFS='_' read -r -a array <<< "${underscored_student_dir}"
		
		shopt -s extdebug
		f "${array[@]}"
		shopt -u extdebug

		student_name_dir=`echo "$*" "${array[@]:3}"`
		
		file_name=`echo $student_name_dir | tr ' ' '_'`

		mv "${student_dir}" "${file_name}"

		cd ${file_name}
		
		for assignment_file in *; do
			mv "${assignment_file}" "${file_name}".pdf
			mv "${file_name}".pdf ..
		done

		cd ..
		
		rmdir ${file_name}
	fi
done
