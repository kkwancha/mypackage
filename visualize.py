import py3Dmol
import math

def showcoord_fromfile(format, file_paths, addlabel=True, showH=False, style='stick'):
    with open(file_paths, 'r') as file:
        string = file.read()
 
    view = py3Dmol.view(width=400, height=300)
    view.addModel(string, format)
    view.setStyle({style: {'colorscheme': 'Jmol'}})
    if addlabel:
        view.addPropertyLabels('index',
                            {'not': {'elem': 'H'}}, 
                            {'fontSize': 10, 
                                'color':'black'})
                                # 'fontColor': 'black',
                                # 'showBackground': False});
        if showH:
            view.addPropertyLabels('index',
                        {}, 
                        {'fontSize': 10, 
                            'color':'black'})
                            # 'fontColor': 'black',
                            # 'showBackground': False});
    view.zoomTo()
    return view.show()

def showcoord_fromstr(format, string , addlabel=True, style='stick'):
    view = py3Dmol.view(width=400, height=300)
    view.addModel(string, format)
    view.setStyle({style: {'colorscheme': 'Jmol'}})
    if addlabel:
        view.addPropertyLabels('index',
                            {'not': {'elem': 'H'}}, 
                            {'fontSize': 10, 
                                'color':'black'})
                                # 'fontColor': 'black',
                                # 'showBackground': False});
    view.zoomTo()
    return view.show()

def showcoords_fromfile(format, file_paths, columns=5, size=150, style='stick'):
    # format = 'xyz', 'mol2'
    columns = 5
    w = size
    h = size
    rows = int(math.ceil(float(len(file_paths)) / columns))
    w = w * columns
    h = h * rows
    
    # Initialize Layout
    view = py3Dmol.view(width=w, height=h, linked=False, viewergrid=(rows, columns))
    
    # Initialize starting positions
    x, y = 0, 0
    
    # Loop through XYZ files
    for xyz_file in file_paths:
        with open(xyz_file, 'r') as f:
            xyz_content = f.read()
        
        view.addModel(xyz_content, format, viewer=(x, y))
        view.zoomTo(viewer=(x, y))
        
        # label_coord = {'x':0, 'y':0, 'z':0}
        # view.addLabel('Hi', {'fontColor': 'black', 'fontSize': 10, 'backgroundOpacity': 0, 'position': coords}, viewer=(x, y))
    
        # Update y and x for grid placement
        if y + 1 < columns:  # Fill in columns
            y += 1
        else:
            x += 1
            y = 0
    
    # view.setStyle({'stick': {'colorscheme': 'Jmol'}})
    view.setStyle({style: {'colorscheme': 'Jmol'}})
    
    view.show()
    
def showcoords_fromstr(format, strings, columns=5, size=150, style='stick'):
    # format = 'xyz', 'mol2'
    columns = 5
    w = size
    h = size
    rows = int(math.ceil(float(len(strings)) / columns))
    w = w * columns
    h = h * rows
    
    # Initialize Layout
    view = py3Dmol.view(width=w, height=h, linked=False, viewergrid=(rows, columns))
    
    # Initialize starting positions
    x, y = 0, 0
    
    # Loop through XYZ files
    for string in strings:
        
        view.addModel(string, format, viewer=(x, y))
        view.zoomTo(viewer=(x, y))
        
        # label_coord = {'x':0, 'y':0, 'z':0}
        # view.addLabel('Hi', {'fontColor': 'black', 'fontSize': 10, 'backgroundOpacity': 0, 'position': coords}, viewer=(x, y))
    
        # Update y and x for grid placement
        if y + 1 < columns:  # Fill in columns
            y += 1
        else:
            x += 1
            y = 0
    
    view.setStyle({style: {'colorscheme': 'Jmol'}})
    view.show()

def showvib_fromxyzstr(xyz_str):
    view = py3Dmol.view(width=400, height=300)
    view.addModel(xyz_str,'xyz',{'vibrate': {'frames':10,'amplitude':1}})
    view.setStyle({'stick': {'colorscheme': 'Jmol'}})
    view.addPropertyLabels('index',
                           {'not': {'elem': 'H'}}, 
                           {'fontSize': 10, 
                            'color':'black'})
                            # 'fontColor': 'black',
                            # 'showBackground': False});
    view.setBackgroundColor('0xeeeeee')
    view.animate({'loop': 'backAndForth'})
    view.zoomTo()
    return view.show()