def update_spreadsheet(spreadsheet_name, sweep):
    import gspread
    import pandas as pd
    import numpy as np

    gc = gspread.service_account()
    sh = gc.open(spreadsheet_name)
    
    worksheet_name = str(sweep.temperature)+'mK_'+sweep.filename
    
    print('Updating worksheet {:s}...'.format(worksheet_name))
    
    # try to open an existing worksheet otherwise creates it
    try:
        ws = sh.worksheet(worksheet_name)
    except:
        print("Selected worksheet does not exist. Creating a new one with the name '{:s}'.".format(worksheet_name))
        ws = sh.add_worksheet(worksheet_name, rows=1000, cols=50)

    all_cells = [c.address for c in ws.get_all_cells()]
    ws.format(all_cells, {"horizontalAlignment": "CENTER", "textFormat": {"fontSize": 8}})
    # bold header
    ws.format('A1:U1', {'textFormat': {'bold': True, "fontSize": 8}})

    # table format 
    ws.format("A2:A1000", {"numberFormat": {"type": "NUMBER" ,"pattern": "0"}})
    ws.format("B2:B1000", {"numberFormat": {"type": "NUMBER", "pattern": "0.000"}})
    ws.format("C2:F1000", {"numberFormat": {"type": "SCIENTIFIC", "pattern": "0.000E+0"}})
    ws.format("G2:L1000", {"numberFormat": {"type": "NUMBER", "pattern": "0"}})
    ws.format("M2:N1000", {"numberFormat": {"type": "NUMBER", "pattern": "0.000"}})
    ws.format("O2:P1000", {"numberFormat": {"type": "SCIENTIFIC", "pattern": "0.000E+0"}})
    ws.format("Q2:Q1000", {"numberFormat": {"type": "NUMBER", "pattern": "0.00"}})

    data = {'channel':      [i for i in range(len(sweep.entry))],
            'target_freq':  [tf['target_freq'] for tf in sweep.entry],
            'Re[a]':        ['-' if e['Re[a]'] is None else 'nan' if np.isnan(e['Re[a]'].n) else e['Re[a]'].n for e in sweep.entry], 
            'Re[a] err':    ['-' if e['Re[a]'] is None else 'nan' if np.isnan(e['Re[a]'].s) else e['Re[a]'].s for e in sweep.entry], 
            'Im[a]':        ['-' if e['Im[a]'] is None else 'nan' if np.isnan(e['Im[a]'].n) else e['Im[a]'].n for e in sweep.entry], 
            'Im[a] err':    ['-' if e['Im[a]'] is None else 'nan' if np.isnan(e['Im[a]'].s) else e['Im[a]'].s for e in sweep.entry], 
            'Q_tot':        ['-' if e['Q_tot'] is None else 'nan' if np.isnan(e['Q_tot'].n) else e['Q_tot'].n for e in sweep.entry], 
            'Q_tot err':    ['-' if e['Q_tot'] is None else 'nan' if np.isnan(e['Q_tot'].s) else e['Q_tot'].s for e in sweep.entry], 
            'Q_i':          ['-' if e['Q_i']   is None else 'nan' if np.isnan(e['Q_i'].n)   else e['Q_i'].n for e in sweep.entry], 
            'Q_i err':      ['-' if e['Q_i']   is None else 'nan' if np.isnan(e['Q_i'].s)   else e['Q_i'].s for e in sweep.entry], 
            'Q_c':          ['-' if e['Q_c']   is None else 'nan' if np.isnan(e['Q_c'].n)   else e['Q_c'].n for e in sweep.entry], 
            'Q_c err':      ['-' if e['Q_c']   is None else 'nan' if np.isnan(e['Q_c'].s)   else e['Q_c'].s for e in sweep.entry], 
            'nu_r':         ['-' if e['nu_r']  is None else 'nan' if np.isnan(e['nu_r'].n)  else e['nu_r'].n for e in sweep.entry], 
            'nu_r err':     ['-' if e['nu_r']  is None else 'nan' if np.isnan(e['nu_r'].s)  else e['nu_r'].s for e in sweep.entry], 
            'phi_0':        ['-' if e['phi_0'] is None else 'nan' if np.isnan(e['phi_0'].n) else e['phi_0'].n for e in sweep.entry], 
            'phi_0 err':    ['-' if e['phi_0'] is None else 'nan' if np.isnan(e['phi_0'].s) else e['phi_0'].s for e in sweep.entry],
            'reduced_chi2': ['-' if e['reduced_chi2'] is None else e['reduced_chi2'] for e in sweep.entry]}

    # dataframe
    df = pd.DataFrame(data=data)

    # write on the gspreadsheet
    ws.update([df.columns.values.tolist()] + df.values.tolist())
    
    print('Updated.')
    
    

"""
this function is deprecated. Use write_column or write_dataframe instead.
"""
def write_electrical_phase_responsivity(spreadsheet_name, index, responsivity):
    import gspread
    import pandas as pd

    gc = gspread.service_account()
    sh = gc.open(spreadsheet_name)
    
    worksheet_name = 'electrical_phase_responsivity'
    
    # try to open an existing worksheet
    try:
        ws = sh.worksheet(worksheet_name)
    except:
        print("Selected worksheet does not exist. Creating a new one with the name '{:s}'.".format(worksheet_name))
        ws = sh.add_worksheet(worksheet_name, rows=1000, cols=50)
    
    # bold header
    ws.format('A1:U1', {'textFormat': {'bold': True, "fontSize": 10}})
    # table format 
    ws.format("A2:A1000", {"numberFormat": {"type": "NUMBER" ,"pattern": "0"}})
    ws.format("B2:C1000", {"numberFormat": {"type": "SCIENTIFIC", "pattern": "0.000E+0"}})
        
    data = {'index':      index,
            'responsivity [rad/W]':  [r.n for r in responsivity],
            'responsivity err [rad/W]':  [r.s for r in responsivity]}
    
    # dataframe
    df = pd.DataFrame(data=data)

    # write on the gspreadsheet
    ws.update([df.columns.values.tolist()] + df.values.tolist())
    
    
    
    
    
def write_column(spreadsheet_name, worksheet_name, header, coordinates, data_column):
    import gspread
    import pandas as pd

    gc = gspread.service_account()
    sh = gc.open(spreadsheet_name)
    
    # try to open an existing worksheet
    try:
        print("Updating '{:s}' -> '{:s}' -> '{:s}'...".format(spreadsheet_name, worksheet_name, header))
        ws = sh.worksheet(worksheet_name)
    except:
        print("Selected worksheet does not exist.")
    
    row, col = gspread.utils.a1_to_rowcol(coordinates)
    
    column_size = len(data_column)
    
    pos1 = gspread.utils.rowcol_to_a1(row+1, col)
    pos2 = gspread.utils.rowcol_to_a1(row+1+column_size, col)
    
    # bold header
    ws.format(coordinates, {'textFormat': {'bold': True, "fontSize": 9}})
    # column format 
    ws.format('{:s}:{:s}'.format(pos1, pos2), {"numberFormat": {"type": "SCIENTIFIC", "pattern": "0.000E+0"}})

    df = pd.DataFrame(data={header: data_column})
    
    # write on the gspreadsheet
    ws.update(coordinates, [df.columns.values.tolist()] + df.values.tolist())
    
    
    
    
def write_dataframe(spreadsheet_name, dataframe, coordinates, worksheet_name=None):
    import gspread

    gc = gspread.service_account()
    sh = gc.open(spreadsheet_name)
    
    # try to open an existing worksheet
    if worksheet_name != None:
        try:
            print("Updating '{:s}' -> '{:s}'...".format(spreadsheet_name, worksheet_name))
            ws = sh.worksheet(worksheet_name)
        except:
            print("Selected worksheet does not exist.")
    else:
        from datetime import datetime
        worksheet_name = datetime.now().strftime("%Y%m%d_%H%M%S")
        print("Updating '{:s}' -> '{:s}'...".format(spreadsheet_name, worksheet_name))
        ws = sh.add_worksheet(worksheet_name, rows=1000, cols=50)
        
    
    row, col = gspread.utils.a1_to_rowcol(coordinates)
    num_rows, num_cols = dataframe.shape
    
    
    # bold header
    if num_cols > 1:
        left_pos = gspread.utils.rowcol_to_a1(row, col)
        right_pos = gspread.utils.rowcol_to_a1(row, col+num_cols)
        ws.format('{:s}:{:s}'.format(left_pos, right_pos), {'textFormat': {'bold': True, "fontSize": 9}})
    else:
        ws.format('{:s}'.format(coordinates), {'textFormat': {'bold': True, "fontSize": 9}})
        
    # column format
    if num_rows > 2:
        upper_left_pos = gspread.utils.rowcol_to_a1(row+1, col)
        lower_right_pos = gspread.utils.rowcol_to_a1(row+1+num_rows, col+num_cols)
        ws.format('{:s}:{:s}'.format(upper_left_pos, lower_right_pos), {"numberFormat": {"type": "SCIENTIFIC", "pattern": "0.000E+0"}})

    # write on the gspreadsheet
    ws.update(coordinates, [dataframe.columns.values.tolist()] + dataframe.values.tolist())
    