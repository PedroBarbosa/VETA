import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns 


def set_style():
    
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")
    
    # This sets reasonable defaults for font size for
    # a figure that will go in a paper
    sns.set_context("paper")
    
    # Set the font to be serif, rather than sans
    sns.set(font='serif')
    
    # Make the background white, and specify the
    # specific font family
    sns.set_style("white", {
        "font.family": "serif",
        "font.serif": ["Times", "Palatino", "serif"]
    })


def set_size(fig, s=10):
    fig.set_size_inches(8, s/3+1)
    plt.tight_layout()


def format_fn(tick_val, tick_pos):
    return abs(tick_val)


def nonna(col):
    return col.isnull().sum() == col.shape[0]
