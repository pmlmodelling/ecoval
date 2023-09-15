import os

import glob
import os
for ff in glob.glob("matched/*.csv"):
    new_ff = "../../" + ff
    variable = os.path.basename(ff).split("_")[1].replace(".csv", "")
    if not os.path.exists(f'book/notebooks/{variable}.ipynb'):
        if variable == "ph":
            Variable = "pH" 
        else:
            Variable = variable.title()
        file1 = 'book/template.ipynb'
        count = 0
        with open(file1, 'r') as file :
          filedata = file.read()
        
        # Replace the target string
        filedata = filedata.replace('template_file_name', new_ff)
        filedata = filedata.replace('template_title', Variable)
        
        # Write the file out again
        with open(f'book/notebooks/{variable}.ipynb', 'w') as file:
          file.write(filedata)
      



os.system("jupyter-book build book/")
