import seaborn as sns

def set_style():
  sns.reset_defaults()
  sns.set_context(
      "notebook", rc={
        "font.size":15,
        "axes.titlesize":20,
        "axes.labelsize":15
      }
  )   
  sns.set_style(
    "darkgrid", 
        {
          'axes.labelcolor': '0',
          'text.color': '0',
          'xtick.color': '0',
          'ytick.color': '0',
          'xtick.bottom': True,
        }
  )
  return True

def get_palette(n, r=False):
    '''returns a suitable set of colors
    n=number of requested colors
    '''
  
    bit_palette = [
      '#817',
      '#a35',
      '#c66',
      '#e94',
      '#ed0',
      '#9d5',
      '#4d8',
      '#2cb',
      '#0bc',
      '#09c',
      '#36b',
      '#639',
    ]
  
    mod=12//n
    pal = [x for n,x in enumerate(bit_palette) if n%mod == 0]
    if r:
        return pal[::-1]
    return pal