from __future__ import annotations

import json

import customtkinter

from . import PrimeEdit, design_guides, design_prime_edit, rank_guides, resolve_nuclease


class App(customtkinter.CTk):
    """Tkinter frontend preserving the original visual hierarchy and typography."""

    def __init__(self):
        super().__init__()
        self.title('gRNA Library')
        self.geometry('1200x1200')
        self.grid_columnconfigure((0, 1), weight=1)
        self.guides = []
        self.current_spec = resolve_nuclease()
        self.reinstate_gRNA_elements()

    def clear_frames(self):
        for name in ('input_frame', 'output_frame'):
            frame = getattr(self, name, None)
            if frame is not None:
                for widget in frame.winfo_children():
                    widget.destroy()

    def _title(self):
        self.title_1 = customtkinter.CTkLabel(self, text='Guide RNA Library', font=('Arial', 120))
        self.title_1.grid(row=0, column=0, pady=12, padx=10, sticky='ew')

    def reinstate_gRNA_elements(self):
        self.clear_frames()
        self._title()
        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.input_frame.grid_columnconfigure((0, 1), weight=1)
        self.entry1 = customtkinter.CTkTextbox(self.input_frame, font=('Helvetica', 30, 'bold'))
        self.entry1.grid(row=0, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.entry1.insert('0.0', "Input gene sequence from 5' to 3' here")
        self.cas_var = customtkinter.IntVar()
        self.cas_box = customtkinter.CTkCheckBox(
            self.input_frame, text='Not Sp Cas9?', variable=self.cas_var, font=('Helvetica', 30, 'bold')
        )
        self.cas_box.grid(row=1, column=0, pady=12, padx=10, sticky='nsew')
        self.cas_var.trace_add('write', self.custom_cas_box)
        self.button = customtkinter.CTkButton(
            self.input_frame, text='Run', command=self.receiver, font=('Helvetica', 50, 'bold')
        )
        self.button.grid(row=3, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
        self.prime_button = customtkinter.CTkButton(
            self.input_frame, text='Prime edit', command=self.reinstate_prime_elements, font=('Helvetica', 30, 'bold')
        )
        self.prime_button.grid(row=4, column=0, pady=12, padx=10, sticky='ew', columnspan=2)

        self.output_frame = customtkinter.CTkFrame(self)
        self.output_frame.grid(row=2, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.output_frame.grid_columnconfigure((0, 1), weight=1)
        self.out_explanation = customtkinter.CTkLabel(
            self.output_frame, text='Output:{gRNA sequence, position, PAM}', font=('Helvetica', 30, 'bold')
        )
        self.out_explanation.grid(row=0, column=0, pady=12, padx=10, sticky='ew')
        self.output1 = customtkinter.CTkTextbox(self.output_frame)
        self.output1.grid(row=1, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.button2 = customtkinter.CTkButton(
            self.output_frame, text='Rank Candidates', command=self.reinstate_ranking_elements, font=('Helvetica', 50, 'bold')
        )
        self.button2.grid(row=2, column=0, pady=12, padx=10, sticky='ew', columnspan=2)

    def custom_cas_box(self, *_):
        if self.cas_var.get() == 1:
            if not hasattr(self, 'custom_cas') or not self.custom_cas.winfo_exists():
                self.custom_cas = customtkinter.CTkEntry(
                    self.input_frame, font=('Helvetica', 30, 'bold'),
                    placeholder_text='Enter the PAM sequence with IUPAC symbols. Ex: NGG'
                )
                self.custom_cas.grid(row=2, column=0, pady=12, padx=10, sticky='nsew', columnspan=2)
        elif hasattr(self, 'custom_cas') and self.custom_cas.winfo_exists():
            self.custom_cas.destroy()

    def _spec(self):
        pam = None
        if self.cas_var.get() == 1 and hasattr(self, 'custom_cas'):
            pam = self.custom_cas.get().strip() or 'NGG'
        return resolve_nuclease('SpCas9', pam=pam)

    def receiver(self):
        seq = self.entry1.get('0.0', 'end').strip()
        self.current_spec = self._spec()
        self.guides = design_guides(seq, self.current_spec)
        lines = [f"{g.spacer} , {g.strand}:{g.start}-{g.end} , {g.pam}" for g in self.guides]
        self.output1.delete('0.0', 'end')
        self.output1.insert('0.0', '\n'.join(lines))

    def reinstate_ranking_elements(self):
        self.clear_frames()
        self._title()
        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.input_frame.grid_columnconfigure((0, 1), weight=1)
        self.entry1 = customtkinter.CTkTextbox(self.input_frame, font=('Helvetica', 30, 'bold'))
        self.entry1.grid(row=0, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.entry1.insert('0.0', "Input the genome sequence to be screened from 5' to 3' here (optional)")
        self.button = customtkinter.CTkButton(
            self.input_frame, text='Run', command=self.rank_candidates, font=('Helvetica', 50, 'bold')
        )
        self.button.grid(row=3, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
        self.output_frame = customtkinter.CTkFrame(self)
        self.output_frame.grid(row=2, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.label_positive = customtkinter.CTkLabel(
            self.output_frame, text='gRNA sequences Ranking and component scores:', font=('Helvetica', 30, 'bold')
        )
        self.label_positive.grid(row=0, column=0, pady=12, padx=10, sticky='nsew')
        self.output1 = customtkinter.CTkTextbox(self.output_frame)
        self.output1.grid(row=1, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.button_3 = customtkinter.CTkButton(
            self.output_frame, text='Find new gRNAs', command=self.reinstate_gRNA_elements, font=('Helvetica', 50, 'bold')
        )
        self.button_3.grid(row=2, column=0, pady=12, padx=10, sticky='ew', columnspan=2)

    def rank_candidates(self):
        genome = self.entry1.get('0.0', 'end').strip()
        if genome.startswith('Input the genome sequence'):
            genome = ''
        ranked = rank_guides(self.guides, spec=self.current_spec, genome_sequence=genome or None)
        lines = [f"#{g.rank} {g.spacer}  {json.dumps(g.scores)}" for g in ranked]
        self.output1.delete('0.0', 'end')
        self.output1.insert('0.0', '\n'.join(lines))

    def reinstate_prime_elements(self):
        self.clear_frames()
        self._title()
        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.input_frame.grid_columnconfigure((0, 1), weight=1)
        self.entry1 = customtkinter.CTkTextbox(self.input_frame, font=('Helvetica', 30, 'bold'))
        self.entry1.grid(row=0, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.entry1.insert('0.0', "Input reference sequence from 5' to 3' here")
        self.position = customtkinter.CTkEntry(self.input_frame, placeholder_text='0-based edit position', font=('Helvetica', 30, 'bold'))
        self.position.grid(row=1, column=0, pady=12, padx=10, sticky='ew')
        self.ref = customtkinter.CTkEntry(self.input_frame, placeholder_text='REF allele', font=('Helvetica', 30, 'bold'))
        self.ref.grid(row=1, column=1, pady=12, padx=10, sticky='ew')
        self.alt = customtkinter.CTkEntry(self.input_frame, placeholder_text='ALT allele', font=('Helvetica', 30, 'bold'))
        self.alt.grid(row=2, column=0, pady=12, padx=10, sticky='ew')
        self.button = customtkinter.CTkButton(self.input_frame, text='Run', command=self.run_prime, font=('Helvetica', 50, 'bold'))
        self.button.grid(row=3, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
        self.output_frame = customtkinter.CTkFrame(self)
        self.output_frame.grid(row=2, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.output1 = customtkinter.CTkTextbox(self.output_frame)
        self.output1.grid(row=0, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.button_3 = customtkinter.CTkButton(self.output_frame, text='Find new gRNAs', command=self.reinstate_gRNA_elements, font=('Helvetica', 50, 'bold'))
        self.button_3.grid(row=1, column=0, pady=12, padx=10, sticky='ew', columnspan=2)

    def run_prime(self):
        seq = self.entry1.get('0.0', 'end').strip()
        candidates = design_prime_edit(seq, PrimeEdit(int(self.position.get()), self.ref.get(), self.alt.get()), resolve_nuclease())
        lines = [
            f"#{c.rank} {c.spacer.spacer} PBS={c.pbs_length} RTT={c.rtt_length} extension={c.extension_sequence} score={c.scores['structural_prime']}"
            for c in candidates[:200]
        ]
        self.output1.delete('0.0', 'end')
        self.output1.insert('0.0', '\n'.join(lines))


def main() -> None:
    app = App()
    app.mainloop()


if __name__ == '__main__':
    main()
