import tkinter as tk
from tkinter import ttk
from tkinter import font
from tkinter import messagebox
from tkinter import filedialog


class Gui():
    def __init__(self):
        self.root = tk.Tk()
        self.style = tk.ttk.Style()
        self.title = "Diffusion calculator"

        self.width = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2
        self.height = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 4

        self.arial = tk.font.Font(family='Times New Roman', size=12, weight=tk.font.BOLD)
        self.info_font = tk.font.Font(family='Times New Roman', size=9, weight=tk.font.BOLD)

        self.note = tk.ttk.Notebook(self.root)
        self.note.pack(fill=tk.BOTH, expand=True)
        self.tabModeling = tk.ttk.Frame(self.note)
        self.tabManual = tk.ttk.Frame(self.note)
        self.tabAbout = tk.ttk.Frame(self.note)
        self.tabModelingInternal = None

        self.enter_experimental_data = tk.BooleanVar()
        self.use_distribution = tk.BooleanVar()
        self.geometry_var = tk.IntVar()
        self.plot_q = tk.BooleanVar()
        self.plot_flux = tk.BooleanVar()
        self.plot_total_flux = tk.BooleanVar()
        self.plot_concentration = tk.BooleanVar()

    def choose_main_options(self, row, column):
        "###-=-=-=-=-=-=-=-=-=- Choose what you want -=-=-=-=-=-=-=-=-=-=-###"
        # does user want to use experimental data
        choose_label = tk.Label(self.tabModelingInternal, text='Please, Choose what you want:', font=self.arial) \
            .grid(row=row, column=column, padx=10, pady=10, sticky=tk.W)
        choose_root = tk.Label(self.tabModelingInternal)
        choose_root.grid(row=row + 1, column=column, padx=10, sticky=tk.W)
        self.enter_experimental_data.set(0)
        tk.Checkbutton(choose_root, text="Enter experimental data",
                       variable=self.enter_experimental_data,
                       onvalue=1, offvalue=0).pack(anchor=tk.W, padx=10)

        # does user want to use size distribution
        self.use_distribution.set(0)
        tk.Checkbutton(choose_root, text="Use distribution",
                       variable=self.use_distribution,
                       onvalue=1, offvalue=0).pack(anchor=tk.W, padx=10)

    def choose_geometry(self, row, column):
        tk.Label(self.tabModelingInternal, text='Please, choose the geometry:', font=self.arial). \
            grid(row=row, column=column, padx=30, pady=10, sticky=tk.W)
        geometry_root = tk.Label(self.tabModelingInternal)
        geometry_root.grid(row=row + 1, column=column, padx=30, sticky=tk.W)
        self.geometry_var.set(0)
        rbutton_slab = tk.Radiobutton(geometry_root,
                                      text='slab',
                                      variable=self.geometry_var,
                                      value=1,
                                      command=lambda: [self.set_geometry_var(1)])
        rbutton_cylinder = tk.Radiobutton(geometry_root,
                                          text='cylinder',
                                          variable=self.geometry_var,
                                          value=2,
                                          command=lambda: [self.set_geometry_var(2)])
        rbutton_sphere = tk.Radiobutton(geometry_root,
                                        text='sphere',
                                        variable=self.geometry_var,
                                        value=3,
                                        command=lambda: [self.set_geometry_var(3)])
        rbutton_slab.grid(row=0, column=0, padx=10, sticky=tk.W)
        rbutton_cylinder.grid(row=0, column=1, padx=10, sticky=tk.W)
        rbutton_sphere.grid(row=0, column=2, padx=10, sticky=tk.W)

    def set_geometry_var(self, val):
        self.geometry_var.set(val)

    def choose_what_to_plot(self, row, column):
        tk.Label(self.tabModelingInternal, text='Please, Choose what to plot:', font=self.arial). \
            grid(row=row, column=column, padx=10, pady=10, sticky=tk.W)
        plot_root = tk.Label(self.tabModelingInternal).grid(row=row + 1, column=column, padx=10, sticky=tk.W)
        self.plot_q.set(0)
        tk.Checkbutton(plot_root, text="total release",
                       variable=self.plot_q,
                       onvalue=1, offvalue=0).grid(row=0, column=0, padx=5, sticky=tk.W)
        self.plot_flux.set(0)
        tk.Checkbutton(plot_root, text="flux",
                       variable=self.plot_flux,
                       onvalue=1, offvalue=0).grid(row=0, column=1, padx=5, sticky=tk.W)
        self.plot_total_flux.set(0)
        tk.Checkbutton(plot_root, text="total flux",
                       variable=self.plot_total_flux,
                       onvalue=1, offvalue=0).grid(row=1, column=0, padx=5, sticky=tk.W)
        self.plot_concentration.set(0)
        tk.Checkbutton(plot_root, text="concentration",
                       variable=self.plot_concentration,
                       onvalue=1, offvalue=0).grid(row=1, column=1, padx=5, sticky=tk.W)

    def help_info(self, text):
        note = tk.Toplevel(self.root)
        xx = (note.winfo_screenwidth() - note.winfo_reqwidth()) / 2
        yy = (note.winfo_screenheight() - note.winfo_reqheight()) / 2
        note.wm_geometry("+%d+%d" % (xx, yy))
        info = tk.Label(note, text=text)
        info.pack()

    def initial_values(self, row, column):
        tk.Label(self.tabModelingInternal, text='Please, enter the initial values.', font=self.arial).\
            grid(row=3, column=0, columnspan=3, pady=5)
        tk.Label(self.tabModelingInternal, text='Enter loaded \n amount of drug (M0):').grid(row=8, column=0)
        ttk.Button(self.tabModelingInternal, text="?", width=3, command=lambda: [self.help_info(
            '''Estimated amount of loaded drug. 
You can use the relative value 
in the persentage (%), then the 
data will be presented on the
graph also in a percentage.''')]).grid(row=8, column=1)
        tk.Label(tab1, text='Enter amount of drug \n on the surface (Q0):').grid(row=9, column=0)
        ttk.Button(tab1, text="?", width=3, command=lambda: [Help_info(
            '''Estimated amount of drug, located
on the vehicle surface.
If you do not know this value,
you can set it to zero (0).
If the M0 was in the %, this value
also must be in the %.''')]).grid(row=9, column=1)
        entry_td = Entry(tab1)
        entry_td.grid(row=7, column=2)
        entry_m0 = Entry(tab1)
        entry_m0.grid(row=8, column=2)
        entry_q0 = Entry(tab1)
        entry_q0.grid(row=9, column=2)

    def tab_modeling(self, row, column):
        self.tabModeling.pack(fill=tk.BOTH, expand=True)
        canvas = tk.Canvas(self.tabModeling)
        vbar = tk.Scrollbar(self.tabModeling, orient=tk.VERTICAL, command=canvas.yview)
        vbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.configure(yscrollcommand=vbar.set, scrollregion=canvas.bbox("all"))
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tabModelingInternal = tk.Frame(canvas)
        canvas.create_window((0, 0), window=self.tabModelingInternal, anchor=tk.NW)
        self.note.add(self.tabModeling, text="Modeling")

        "###-=-=-=-=-=-=-=-=-=- Choose what you want -=-=-=-=-=-=-=-=-=-=-###"
        self.choose_main_options(0, 0)
        "###-=-=-=-=-=-=- Choose the geometry -=-=-=-=-=-=-###"
        self.choose_geometry(2, 1)
        "###-=-=-=-=-=-=- Choose what to plot -=-=-=-=-=-=-###"
        self.choose_what_to_plot(2, 0)


    def tab_manual(self, row, column):
        self.tabManual.pack(fill=tk.BOTH, expand=True)
        self.note.add(self.tabManual, text="Manual")

        tk.Label(self.tabManual, text='Step-by-step Instructions', font=self.arial). \
            grid(row=0, column=0, padx=10, sticky=tk.N, pady=10)
        tk.Label(self.tabManual, text='(1) Under "choose what you want" section, decide whether you have ' +
                                      'experimental data, or a known distribution and if you want to plot ' +
                                      'graphs by changing the truth value of the appropriate variable.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual, text='(2) Under the "choose geometry" section, enter the geometry of the ' +
                                      'problem, which results you want to plot and some style choices for the plots.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 1, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual, text='(3) Under the "define parameters" section, enter the units of your input ' +
                                      'parameters and the units that you want to be graphed (output), then ' +
                                      'enter all parameters.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 2, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual, text='(4) If you have experimental data. Under "Enter experimental data" ' +
                                      'section, import your data from an excel document (up to 5 files), it ' +
                                      'will be converted into a numpy array with transposed rows/columns. ' +
                                      'Then assign the data to either flux, total flux, concentration or ' +
                                      'release, (up to 5 of each) as appropriate.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 3, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual, text='Notes About Parameter Choices', font=self.arial). \
            grid(row=row + 4, column=column, padx=10, sticky=tk.N, pady=10)
        tk.Label(self.tabManual, text='- Don\'t make mass and mass_bound too small.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 5, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual, text='- The tolerance should be much smaller for slabs than other geometries.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 6, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual, text='- The larger the value chosen for N is, ti must chosen sufficiently large and ' +
                                      'tolerance sufficiently low.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 7, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual,
                 text='- Make sure SD is appropriately small for the choice of R0, if SD is too large ' +
                      'there will be particles with negative radius, the equality on line 287 ' +
                      'will always hold and the solver will not run, and an error will be ' +
                      'produced on line 395.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 8, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual, text='- When v is made too large, discontinuities in the flux can become bad ' +
                                      'v should not be more than 1 order of magnitude larger than D.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 9, column=column, padx=10, sticky=tk.W, pady=2)
        tk.Label(self.tabManual,
                 text='- Diffusion coefficients are much smaller in nanoparticles/fibres than in slabs.',
                 wraplength=self.root.winfo_reqwidth() * 5). \
            grid(row=row + 10, column=column, padx=10, sticky=tk.W, pady=2)

    def tab_about(self):
        self.tabAbout.pack(fill=tk.BOTH, expand=True)
        self.note.add(self.tabAbout, text="About")

        tk.Label(self.tabAbout, text="This program based on the compartmental and size distribution models " +
                                     "and was developed by the Griffith University and " +
                                     "Tomsk Polytechnic University research team \n " +
                                     "(Anissimov Y.G., Spiridonova T.I., Marriott R. and Petlin D.G.)",
                 wraplength=self.root.winfo_reqwidth() * 2).pack(pady=10)
        tk.Label(self.tabAbout,
                 text='The contact e-mails:\ny.anissimov@griffith.edu.au\nrory.marriott@griffithuni.edu.au',
                 wraplength=self.root.winfo_reqwidth() * 2, anchor='w').pack(pady=10)

    def mainloop(self, command):
        self.tab_modeling()
        self.tab_manual(1, 0)
        self.tab_about()

        self.style.configure('Calculate.TButton', background='#3CB371', font=self.arial, width=20)
        calculate_button = ttk.Button(self.tabModelingInternal,
                                      text='Calculate',
                                      command=command,
                                      style='Calculate.TButton')
        calculate_button.grid(row=21, column=2, sticky=tk.E)
        emply_label = tk.Label(self.tabModelingInternal, text='')
        emply_label.grid(row=22, column=0)

        self.root.mainloop()
