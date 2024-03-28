file_names = ['figure1.py' , 'figure2.py' , 'figure3.py' ,
              'figure4.py' , 'figure5.py' , 'figure6.py' ,
              'figure7.py' , 'figure8.py' , 'figure9.py' ,
              'figureA1.py','figureA2.py' , 'figureA3.py',
              'figureB1_alpha.py', 'figureB1_slope.py' ]

figures_directory = './Figures_py/'

for file_name in file_names:
    with open(figures_directory + file_name, 'r') as file:
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