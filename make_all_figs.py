import warnings
warnings.filterwarnings("ignore")

file_names = ['Figure2.py' , 'Figure3.py' , 'Figure4.py',
              'appendix_A1.py', 'appendix_B.py' ]

for file_name in file_names:
    with open(file_name, 'r') as file:
        print('\n')
        to_output = f'#### Starting: {file_name} ####'
        print('#'*len(to_output))
        print(to_output)
        print('#'*len(to_output))
        print('\n')
        exec(file.read())
        print('\n')
        print('!'*len(to_output))
        print(f'!!!! Finished: {file_name} !!!!')
        print('!'*len(to_output))
