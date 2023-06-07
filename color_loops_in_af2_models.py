import pandas as pd

# read csv file output of get_loop_residues.py
df = pd.read_csv('diversity_panel_loops.csv')

# remove unnecessary columns
df = df[['Seq', 'loop_1_start', 'loop_1_end', \
        'loop_2_start', 'loop_2_end', \
        'loop_3_start', 'loop_3_end', \
        'loop_4_start', 'loop_4_end', \
        'loop_5_start', 'loop_5_end', \
        'loop_6_start', 'loop_6_end', \
        ]]
df

# create pymol selections
starts = ['loop_1_start', 'loop_2_start', 'loop_3_start', 'loop_4_start', 'loop_5_start', 'loop_6_start']
ends = ['loop_1_end', 'loop_2_end', 'loop_3_end', 'loop_4_end', 'loop_5_end', 'loop_6_end']
colors = ['red', 'orange', 'yellow', 'blue', 'purple', 'gray']

with open('loop_selections.pml', 'w') as pml:
    for i in range(len(df)):        
        for start, end, color in zip(starts, ends, colors):
            # loop 6 is on the other chain
            if (start != 'loop_6_start'):
                chain = 'chain A'
            else:
                chain = 'chain B'
            
            selection = 'select ' + df.loc[i, 'Seq'] + '_' + start[0:6] + ', '  \
                + df.loc[i, 'Seq'] + ' and ' + chain + ' and resi ' + str(df.loc[i, start]) \
                + '-' + str(df.loc[i, end])
            
            color = 'color ' + color + ', ' + df.loc[i, 'Seq'] + '_' + start[0:6]
            #print(selection)
            #print(color + '\n')
            pml.write(selection + '\n')
            pml.write(color + '\n')
