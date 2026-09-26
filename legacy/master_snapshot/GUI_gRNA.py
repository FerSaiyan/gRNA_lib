import customtkinter
import gRNA_lib_GUI as glib


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.title("gRNA Library")
        self.geometry('1200x1200')
        self.grid_columnconfigure((0, 1), weight=1)

        self.reinstate_gRNA_elements()

    def reinstate_gRNA_elements(self):
        self.clear_frames()
        self.title_1 = customtkinter.CTkLabel(self, text='Guide RNA Library', font=('Arial', 120))
        self.title_1.grid(row=0, column=0, pady=12, padx=10, sticky='ew')

        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.input_frame.grid_columnconfigure((0, 1), weight=1)

        self.entry1 = customtkinter.CTkTextbox(self.input_frame, font=('Helvetica', 30, 'bold'))
        self.entry1.grid(row=0, column=0, pady=12, padx=10, sticky= "ew", columnspan=3)
        self.entry1.insert("0.0", "Input gene sequence from 5' to 3' here")

        self.cas_var = customtkinter.IntVar()
        self.cas_box = customtkinter.CTkCheckBox(self.input_frame, text='Not Sp Cas9?', variable=self.cas_var, font=('Helvetica', 30, 'bold'))
        self.cas_box.grid(row=1, column=0, pady=12, padx=10, sticky='nsew')
        self.cas_var.trace('w', self.custom_cas_box)
    
        self.button = customtkinter.CTkButton(self.input_frame, text='Run', command=self.receiver, font=('Helvetica', 50, 'bold'))
        self.button.grid(row=3, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
        
        self.output_frame = customtkinter.CTkFrame(self)
        self.output_frame.grid(row=2, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.output_frame.grid_columnconfigure((0, 1), weight=1)
        
        self.out_explanation = customtkinter.CTkLabel(self.output_frame, text='Output:{gRNA sequence, position, PAM}', font=('Helvetica', 30, 'bold'))
        self.out_explanation.grid(row=0, column=0, pady=12, padx=10, sticky='ew')
        
        self.label_positive = customtkinter.CTkLabel(self.output_frame, text='On the Positive strand:', font=('Helvetica', 30, 'bold'))
        self.label_positive.grid(row=1, column=0, pady=12, padx=10, sticky='nsew')
        self.output1 = customtkinter.CTkTextbox(self.output_frame)
        self.output1.grid(row=2, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        
        self.label_negative = customtkinter.CTkLabel(self.output_frame, text='On the Negative strand:', font=('Helvetica', 30, 'bold'))
        self.label_negative.grid(row=3, column=0, pady=12, padx=10, sticky='nsew')
        self.output2 = customtkinter.CTkTextbox(self.output_frame)
        self.output2.grid(row=4, column=0, pady=12, padx=10, sticky='ew', columnspan=3)

        self.button2 = customtkinter.CTkButton(self.output_frame, text='Rank Candidates', command=self.reinstate_ranking_elements, font=('Helvetica', 50, 'bold')) 
        self.button2.grid(row=5, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
    def clear_frames(self):
        if hasattr(self, 'input_frame'):
            #self.input_frame.destroy()
            for widget in self.input_frame.winfo_children():
                widget.destroy()
        if hasattr(self, 'output_frame'):
            #self.output_frame.destroy()
            for widget in self.output_frame.winfo_children():
                widget.destroy()

    def receiver(self):
        self.output1.delete("0.0", "end")
        self.output2.delete("0.0", "end")
        if hasattr(self, 'custom_cas'):
            globals()["gRNA_list_plus"], globals()["gRNA_list_minus"] = glib.main(self.entry1.get('0.0', 'end'), PAM=self.custom_cas.get())
        else:
            globals()["gRNA_list_plus"], globals()["gRNA_list_minus"] = glib.main(self.entry1.get('0.0', 'end'))
        self.output1.insert("0.0", gRNA_list_plus)
        self.output2.insert("0.0", gRNA_list_minus)

    def custom_cas_box(self, *args):
        if self.cas_var.get() == 1:
            self.custom_cas = customtkinter.CTkEntry(self.input_frame, font=('Helvetica', 30, 'bold'), placeholder_text='Enter the PAM sequence with N for any base. Ex: NGG')
            self.custom_cas.grid(row=2, column=0, pady=12, padx=10, sticky='nsew', columnspan=2)
            #self.custom_cas.insert("0.0", 'Enter the PAM sequence with N for any base. Ex: NGG')
            #self.custom_cas.grid_columnconfigure((0, 1), weight=1)
        else:
            if hasattr(self, 'custom_cas'):
                self.custom_cas.destroy()

    def reinstate_ranking_elements(self):
        self.clear_frames()
        self.title_1 = customtkinter.CTkLabel(self, text='Guide RNA Library', font=('Arial', 120))
        self.title_1.grid(row=0, column=0, pady=12, padx=10, sticky='ew')

        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.input_frame.grid_columnconfigure((0, 1), weight=1)

        self.entry1 = customtkinter.CTkTextbox(self.input_frame, font=('Helvetica', 30, 'bold'))
        self.entry1.grid(row=0, column=0, pady=12, padx=10, sticky= "ew", columnspan=3)
        self.entry1.insert("0.0", "Input the genome sequence to be screened from 5' to 3' here")
    
        self.button = customtkinter.CTkButton(self.input_frame, text='Run', command=self.rank_candidates, font=('Helvetica', 50, 'bold'))
        self.button.grid(row=3, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
        
        self.output_frame = customtkinter.CTkFrame(self)
        self.output_frame.grid(row=2, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.output_frame.grid_columnconfigure((0, 1), weight=1)
        
        self.label_positive = customtkinter.CTkLabel(self.output_frame, text='gRNA sequences Ranking and Overall Score:', font=('Helvetica', 30, 'bold'))
        self.label_positive.grid(row=0, column=0, pady=12, padx=10, sticky='nsew')
        self.out_explanation = customtkinter.CTkLabel(self.output_frame, text='Output:{gRNA sequence, Number of miss-matches, Overall score}', font=('Helvetica', 30, 'bold'))
        self.out_explanation.grid(row=1, column=0, pady=12, padx=10, sticky='ew')
        self.output1 = customtkinter.CTkTextbox(self.output_frame)
        self.output1.grid(row=2, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        
        self.button_3 = customtkinter.CTkButton(self.output_frame, text='Find new gRNAs', command=self.reinstate_gRNA_elements, font=('Helvetica', 50, 'bold'))
        self.button_3.grid(row=3, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
    
    def rank_candidates(self):
        ranking = glib.gRNA_ranking((globals()["gRNA_list_plus"]+globals()["gRNA_list_minus"]), self.entry1.get('0.0', 'end'), True if self.cas_var.get() == 1 else False)
        self.output1.delete("0.0", "end")
        aux_rank_list = [(str_x[0] + ' , ' + str(str_x[1])) for str_x in ranking]
        self.output1.insert("0.0", '\n'.join(aux_rank_list))

app = App()
app.mainloop()