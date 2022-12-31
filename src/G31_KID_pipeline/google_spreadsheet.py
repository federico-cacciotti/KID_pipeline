def update_spreadsheet(spreadsheet_name, sweep):
    import gspread
    import pandas as pd
    import numpy as np

    gc = gspread.service_account()
    sh = gc.open(spreadsheet_name)
    
    worksheet_name = str(sweep.temperature)+'mK_'+sweep.filename
    
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