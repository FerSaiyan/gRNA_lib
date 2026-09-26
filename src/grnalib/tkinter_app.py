from __future__ import annotations

import json
import threading
from tkinter import filedialog

import customtkinter

from . import PrimeEdit, design_guides, design_prime_edit, rank_guides, resolve_nuclease
from .genome_sources import (
    download_ncbi_genome,
    fasta_records,
    fetch_ensembl_region,
    list_ncbi_sequences,
    read_fasta_record,
    search_ncbi_assemblies,
)
from .reference_workflows import build_or_reuse_ncbi_index


class App(customtkinter.CTk):
    """Tkinter frontend preserving the original visual hierarchy and typography."""

    def __init__(self):
        super().__init__()
        self.title('gRNA Library')
        self.geometry('1200x1200')
        self.grid_columnconfigure((0, 1), weight=1)
        self.guides = []
        self.current_spec = resolve_nuclease()
        self.cached_ncbi = None
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

    def _status(self, text: str) -> None:
        label = getattr(self, 'source_status', None)
        if label is not None and label.winfo_exists():
            label.configure(text=text)

    def _put_sequence(self, sequence: str, label: str) -> None:
        self.entry1.delete('0.0', 'end')
        self.entry1.insert('0.0', sequence)
        self._status(f'Loaded {label} ({len(sequence):,} bp)')

    def _choose_record(self, names: list[str], title: str = 'Choose FASTA record') -> str | None:
        if len(names) == 1:
            return names[0]
        preview = ', '.join(names[:12])
        if len(names) > 12:
            preview += ', ...'
        dialog = customtkinter.CTkInputDialog(
            text=f'This FASTA has {len(names)} records. Enter one record name:\n{preview}',
            title=title,
        )
        choice = (dialog.get_input() or '').strip()
        return choice or None

    def choose_fasta(self):
        path = filedialog.askopenfilename(
            title='Choose FASTA file',
            filetypes=[
                ('FASTA files', '*.fa *.fasta *.fna *.fas'),
                ('Text files', '*.txt'),
                ('All files', '*.*'),
            ],
        )
        if not path:
            return
        try:
            infos = fasta_records(path)
            record = self._choose_record([info.name for info in infos])
            if record is None:
                self._status('FASTA selection cancelled')
                return
            _, sequence = read_fasta_record(path, record)
            self._put_sequence(sequence, f'{record} from {path}')
        except Exception as exc:
            self._status(str(exc))

    def _add_sequence_sources(
        self,
        start_row: int,
        *,
        include_ncbi: bool = True,
        include_ensembl: bool = True,
    ) -> int:
        row = start_row
        self.fasta_button = customtkinter.CTkButton(
            self.input_frame,
            text='Choose FASTA file',
            command=self.choose_fasta,
            font=('Helvetica', 24, 'bold'),
        )
        self.fasta_button.grid(row=row, column=0, pady=8, padx=10, sticky='ew', columnspan=2)
        row += 1

        if include_ncbi:
            self.ncbi_organism = customtkinter.CTkEntry(
                self.input_frame,
                placeholder_text='Organism, e.g. Homo sapiens',
                font=('Helvetica', 22, 'bold'),
            )
            self.ncbi_organism.grid(row=row, column=0, pady=8, padx=10, sticky='ew')
            self.ncbi_search_button = customtkinter.CTkButton(
                self.input_frame,
                text='Find NCBI assemblies',
                command=self.search_ncbi_organism,
                font=('Helvetica', 22, 'bold'),
            )
            self.ncbi_search_button.grid(row=row, column=1, pady=8, padx=10, sticky='ew')
            row += 1

            self.ncbi_assembly_map = {}
            self.ncbi_assembly_menu = customtkinter.CTkOptionMenu(
                self.input_frame,
                values=['Search assemblies first'],
                command=lambda _: self.load_ncbi_sequences(),
                state='disabled',
                font=('Helvetica', 20, 'bold'),
            )
            self.ncbi_assembly_menu.grid(row=row, column=0, pady=8, padx=10, sticky='ew', columnspan=2)
            row += 1

            self.ncbi_sequence_map = {'Whole assembly': None}
            self.ncbi_sequence_menu = customtkinter.CTkOptionMenu(
                self.input_frame,
                values=['Whole assembly'],
                state='disabled',
                font=('Helvetica', 20, 'bold'),
            )
            self.ncbi_sequence_menu.grid(row=row, column=0, pady=8, padx=10, sticky='ew', columnspan=2)
            row += 1

            self.ncbi_download_button = customtkinter.CTkButton(
                self.input_frame,
                text='Download / use cached reference',
                command=self.download_selected_ncbi,
                state='disabled',
                font=('Helvetica', 22, 'bold'),
            )
            self.ncbi_download_button.grid(row=row, column=0, pady=8, padx=10, sticky='ew')
            self.ncbi_index_button = customtkinter.CTkButton(
                self.input_frame,
                text='Build / reuse off-target index',
                command=self.build_selected_ncbi_index,
                state='disabled',
                font=('Helvetica', 22, 'bold'),
            )
            self.ncbi_index_button.grid(row=row, column=1, pady=8, padx=10, sticky='ew')
            row += 1

            self.index_pam = customtkinter.CTkEntry(
                self.input_frame,
                placeholder_text='Index PAM',
                font=('Helvetica', 20, 'bold'),
            )
            self.index_pam.insert(0, 'NGG')
            self.index_pam.grid(row=row, column=0, pady=8, padx=10, sticky='ew')
            self.index_spacer = customtkinter.CTkEntry(
                self.input_frame,
                placeholder_text='Spacer length',
                font=('Helvetica', 20, 'bold'),
            )
            self.index_spacer.insert(0, '20')
            self.index_spacer.grid(row=row, column=1, pady=8, padx=10, sticky='ew')
            row += 1
            self.index_pam_side = customtkinter.CTkOptionMenu(
                self.input_frame,
                values=['3prime', '5prime'],
                font=('Helvetica', 20, 'bold'),
            )
            self.index_pam_side.set('3prime')
            self.index_pam_side.grid(row=row, column=0, pady=8, padx=10, sticky='ew', columnspan=2)
            row += 1

            self.ncbi_accession = customtkinter.CTkEntry(
                self.input_frame,
                placeholder_text='Or enter NCBI assembly accession directly',
                font=('Helvetica', 20, 'bold'),
            )
            self.ncbi_accession.grid(row=row, column=0, pady=8, padx=10, sticky='ew')
            self.ncbi_chromosome = customtkinter.CTkEntry(
                self.input_frame,
                placeholder_text='Chromosome (optional)',
                font=('Helvetica', 20, 'bold'),
            )
            self.ncbi_chromosome.grid(row=row, column=1, pady=8, padx=10, sticky='ew')
            row += 1
            self.ncbi_fetch_button = customtkinter.CTkButton(
                self.input_frame,
                text='Fetch accession directly',
                command=self.fetch_ncbi,
                font=('Helvetica', 20, 'bold'),
            )
            self.ncbi_fetch_button.grid(row=row, column=0, pady=8, padx=10, sticky='ew')
            self.ncbi_load_button = customtkinter.CTkButton(
                self.input_frame,
                text='Load cached NCBI record',
                command=self.load_cached_ncbi_record,
                state='disabled',
                font=('Helvetica', 20, 'bold'),
            )
            self.ncbi_load_button.grid(row=row, column=1, pady=8, padx=10, sticky='ew')
            row += 1

        if include_ensembl:
            self.ensembl_species = customtkinter.CTkEntry(
                self.input_frame,
                placeholder_text='Ensembl species, e.g. homo_sapiens',
                font=('Helvetica', 22, 'bold'),
            )
            self.ensembl_species.grid(row=row, column=0, pady=8, padx=10, sticky='ew')
            self.ensembl_region = customtkinter.CTkEntry(
                self.input_frame,
                placeholder_text='Region, e.g. 17:7668402..7687550:1',
                font=('Helvetica', 22, 'bold'),
            )
            self.ensembl_region.grid(row=row, column=1, pady=8, padx=10, sticky='ew')
            row += 1
            self.ensembl_fetch_button = customtkinter.CTkButton(
                self.input_frame,
                text='Fetch Ensembl region',
                command=self.fetch_ensembl,
                font=('Helvetica', 24, 'bold'),
            )
            self.ensembl_fetch_button.grid(row=row, column=0, pady=8, padx=10, sticky='ew', columnspan=2)
            row += 1

        self.source_status = customtkinter.CTkLabel(
            self.input_frame,
            text='Paste DNA, choose FASTA, or fetch a public reference.',
            font=('Helvetica', 16),
            wraplength=1080,
        )
        self.source_status.grid(row=row, column=0, pady=6, padx=10, sticky='ew', columnspan=2)
        return row + 1

    def search_ncbi_organism(self):
        taxon = self.ncbi_organism.get().strip()
        if not taxon:
            self._status('Enter an organism name or NCBI TaxID')
            return
        self._status('Searching NCBI assemblies...')
        self.ncbi_search_button.configure(state='disabled')

        def worker():
            try:
                assemblies = search_ncbi_assemblies(taxon, limit=30)
                self.after(0, lambda assemblies=assemblies: self._finish_ncbi_assemblies(assemblies))
            except Exception as exc:
                button = self.ncbi_search_button
                self.after(0, lambda exc=exc, button=button: self._finish_remote_error(exc, button))

        threading.Thread(target=worker, daemon=True).start()

    def _finish_ncbi_assemblies(self, assemblies):
        if hasattr(self, 'ncbi_search_button') and self.ncbi_search_button.winfo_exists():
            self.ncbi_search_button.configure(state='normal')
        if not assemblies:
            self._status('No current RefSeq assemblies were found for that organism')
            return
        if not hasattr(self, 'ncbi_assembly_menu') or not self.ncbi_assembly_menu.winfo_exists():
            return
        self.ncbi_assembly_map = {}
        labels = []
        for assembly in assemblies:
            tag = ' [reference]' if assembly.is_reference else ' [representative]' if assembly.is_representative else ''
            label = f'{assembly.assembly_name or assembly.accession} — {assembly.accession}{tag}'
            if assembly.synonym:
                label += f' — {assembly.synonym}'
            labels.append(label)
            self.ncbi_assembly_map[label] = assembly
        self.ncbi_assembly_menu.configure(values=labels, state='normal')
        self.ncbi_assembly_menu.set(labels[0])
        self.load_ncbi_sequences()

    def _selected_ncbi(self):
        if not hasattr(self, 'ncbi_assembly_menu'):
            return None, None
        assembly = self.ncbi_assembly_map.get(self.ncbi_assembly_menu.get())
        chromosome = self.ncbi_sequence_map.get(self.ncbi_sequence_menu.get())
        return assembly, chromosome

    def load_ncbi_sequences(self):
        assembly = self.ncbi_assembly_map.get(self.ncbi_assembly_menu.get())
        if assembly is None:
            return
        self._status(f'Loading chromosomes for {assembly.accession}...')

        def worker():
            try:
                records = list_ncbi_sequences(assembly.accession, chromosomes_only=True)
                self.after(0, lambda records=records: self._finish_ncbi_sequences(records))
            except Exception as exc:
                self.after(0, lambda exc=exc: self._status(str(exc)))

        threading.Thread(target=worker, daemon=True).start()

    def _finish_ncbi_sequences(self, records):
        if not hasattr(self, 'ncbi_sequence_menu') or not self.ncbi_sequence_menu.winfo_exists():
            return
        self.ncbi_sequence_map = {'Whole assembly': None}
        labels = ['Whole assembly']
        for record in records:
            accession = record.refseq_accession or record.genbank_accession or record.sequence_name
            label = f'chr{record.chromosome} — {accession} — {record.length:,} bp'
            labels.append(label)
            self.ncbi_sequence_map[label] = record.chromosome
        self.ncbi_sequence_menu.configure(values=labels, state='normal')
        self.ncbi_sequence_menu.set(labels[0])
        self.ncbi_download_button.configure(state='normal')
        self.ncbi_index_button.configure(state='normal')
        self._status('Choose the whole assembly or one chromosome. Downloads are cached locally.')

    def download_selected_ncbi(self):
        assembly, chromosome = self._selected_ncbi()
        if assembly is None:
            self._status('Choose an NCBI assembly first')
            return
        self._status('Downloading / checking NCBI cache...')
        self.ncbi_download_button.configure(state='disabled')

        def worker():
            try:
                result = download_ncbi_genome(
                    assembly.accession,
                    chromosomes=[chromosome] if chromosome else None,
                )
                self.after(0, lambda result=result: self._finish_selected_ncbi_download(result))
            except Exception as exc:
                button = self.ncbi_download_button
                self.after(0, lambda exc=exc, button=button: self._finish_remote_error(exc, button))

        threading.Thread(target=worker, daemon=True).start()

    def _finish_selected_ncbi_download(self, result):
        self.cached_ncbi = result
        if hasattr(self, 'ncbi_download_button') and self.ncbi_download_button.winfo_exists():
            self.ncbi_download_button.configure(state='normal')
        if hasattr(self, 'ncbi_load_button') and self.ncbi_load_button.winfo_exists():
            self.ncbi_load_button.configure(state='normal')
        if len(result.records) == 1 and result.records[0].length <= 10_000_000:
            _, sequence = read_fasta_record(result.fasta_path, result.records[0].name)
            self._put_sequence(sequence, f'{result.records[0].name} from NCBI')
        else:
            state = 'Using cached' if result.cached else 'Downloaded'
            self._status(f'{state} reference at {result.fasta_path}. Large records stay on disk for indexed analysis.')

    def build_selected_ncbi_index(self):
        assembly, chromosome = self._selected_ncbi()
        if assembly is None:
            self._status('Choose an NCBI assembly first')
            return
        pam = self.index_pam.get().strip() or 'NGG'
        try:
            spacer_length = int(self.index_spacer.get())
        except ValueError:
            self._status('Spacer length must be an integer')
            return
        pam_side = self.index_pam_side.get()
        self._status('Building or checking the off-target index...')
        self.ncbi_index_button.configure(state='disabled')

        def worker():
            try:
                result = build_or_reuse_ncbi_index(
                    assembly.accession,
                    chromosomes=[chromosome] if chromosome else None,
                    pam=pam,
                    spacer_length=spacer_length,
                    pam_side=pam_side,
                )
                self.after(0, lambda result=result: self._finish_ncbi_index(result))
            except Exception as exc:
                button = self.ncbi_index_button
                self.after(0, lambda exc=exc, button=button: self._finish_remote_error(exc, button))

        threading.Thread(target=worker, daemon=True).start()

    def _finish_ncbi_index(self, result):
        if hasattr(self, 'ncbi_index_button') and self.ncbi_index_button.winfo_exists():
            self.ncbi_index_button.configure(state='normal')
        state = 'Reused cached' if result.cached else 'Built'
        self._status(f'{state} off-target index: {result.index_prefix}')

    def fetch_ncbi(self):
        accession = self.ncbi_accession.get().strip()
        chromosome = self.ncbi_chromosome.get().strip()
        if not accession:
            self._status('Enter an NCBI GCF_/GCA_ assembly accession')
            return
        self._status('Fetching from NCBI Datasets...')
        self.ncbi_fetch_button.configure(state='disabled')

        def worker():
            try:
                result = download_ncbi_genome(
                    accession,
                    chromosomes=[chromosome] if chromosome else None,
                )
                self.after(0, lambda result=result: self._finish_ncbi(result))
            except Exception as exc:
                button = self.ncbi_fetch_button
                self.after(0, lambda exc=exc, button=button: self._finish_remote_error(exc, button))

        threading.Thread(target=worker, daemon=True).start()

    def _finish_ncbi(self, result):
        self.cached_ncbi = result
        if hasattr(self, 'ncbi_fetch_button') and self.ncbi_fetch_button.winfo_exists():
            self.ncbi_fetch_button.configure(state='normal')
        if hasattr(self, 'ncbi_load_button') and self.ncbi_load_button.winfo_exists():
            self.ncbi_load_button.configure(state='normal')
        if len(result.records) == 1 and result.records[0].length <= 10_000_000:
            _, sequence = read_fasta_record(result.fasta_path, result.records[0].name)
            self._put_sequence(sequence, f'{result.records[0].name} from NCBI')
            return
        record_summary = ', '.join(
            f'{record.name} ({record.length:,} bp)' for record in result.records[:6]
        )
        if len(result.records) > 6:
            record_summary += ', ...'
        self._status(
            f'Downloaded to {result.fasta_path}. Records: {record_summary}. '
            'Use “Load cached NCBI record” to place one in the sequence box.'
        )

    def load_cached_ncbi_record(self):
        result = self.cached_ncbi
        if result is None:
            self._status('No NCBI download is cached in this session')
            return
        try:
            record = self._choose_record([info.name for info in result.records], 'Choose NCBI record')
            if record is None:
                return
            _, sequence = read_fasta_record(result.fasta_path, record)
            self._put_sequence(sequence, f'{record} from NCBI cache')
        except Exception as exc:
            self._status(str(exc))

    def fetch_ensembl(self):
        species = self.ensembl_species.get().strip()
        region = self.ensembl_region.get().strip()
        if not species or not region:
            self._status('Enter both an Ensembl species and genomic region')
            return
        self._status('Fetching from Ensembl...')
        self.ensembl_fetch_button.configure(state='disabled')

        def worker():
            try:
                sequence = fetch_ensembl_region(species, region)
                self.after(
                    0,
                    lambda sequence=sequence, species=species, region=region: self._finish_ensembl(sequence, species, region),
                )
            except Exception as exc:
                button = self.ensembl_fetch_button
                self.after(0, lambda exc=exc, button=button: self._finish_remote_error(exc, button))

        threading.Thread(target=worker, daemon=True).start()

    def _finish_ensembl(self, sequence: str, species: str, region: str):
        if hasattr(self, 'ensembl_fetch_button') and self.ensembl_fetch_button.winfo_exists():
            self.ensembl_fetch_button.configure(state='normal')
        self._put_sequence(sequence, f'{species} {region} from Ensembl')

    def _finish_remote_error(self, exc: Exception, button):
        if button is not None and button.winfo_exists():
            button.configure(state='normal')
        self._status(str(exc))

    def reinstate_gRNA_elements(self):
        self.clear_frames()
        self._title()
        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.input_frame.grid_columnconfigure((0, 1), weight=1)
        self.entry1 = customtkinter.CTkTextbox(self.input_frame, font=('Helvetica', 30, 'bold'))
        self.entry1.grid(row=0, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.entry1.insert('0.0', "Input gene sequence from 5' to 3' here")

        row = self._add_sequence_sources(1, include_ncbi=True, include_ensembl=True)

        self.cas_var = customtkinter.IntVar()
        self.cas_box = customtkinter.CTkCheckBox(
            self.input_frame, text='Not Sp Cas9?', variable=self.cas_var, font=('Helvetica', 30, 'bold')
        )
        self.cas_box.grid(row=row, column=0, pady=12, padx=10, sticky='nsew')
        self.custom_cas_row = row + 1
        self.cas_var.trace_add('write', self.custom_cas_box)
        self.button = customtkinter.CTkButton(
            self.input_frame, text='Run', command=self.receiver, font=('Helvetica', 50, 'bold')
        )
        self.button.grid(row=row + 2, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
        self.prime_button = customtkinter.CTkButton(
            self.input_frame, text='Prime edit', command=self.reinstate_prime_elements, font=('Helvetica', 30, 'bold')
        )
        self.prime_button.grid(row=row + 3, column=0, pady=12, padx=10, sticky='ew', columnspan=2)

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
                self.custom_cas.grid(
                    row=self.custom_cas_row,
                    column=0,
                    pady=12,
                    padx=10,
                    sticky='nsew',
                    columnspan=2,
                )
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

        row = self._add_sequence_sources(1, include_ncbi=True, include_ensembl=False)
        self.source_status.configure(
            text='For short local references, paste DNA or choose FASTA. '
            'Whole-genome work should use an indexed backend.'
        )
        self.button = customtkinter.CTkButton(
            self.input_frame, text='Run', command=self.rank_candidates, font=('Helvetica', 50, 'bold')
        )
        self.button.grid(row=row, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
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

        row = self._add_sequence_sources(1, include_ncbi=False, include_ensembl=True)
        self.position = customtkinter.CTkEntry(
            self.input_frame, placeholder_text='0-based edit position', font=('Helvetica', 30, 'bold')
        )
        self.position.grid(row=row, column=0, pady=12, padx=10, sticky='ew')
        self.ref = customtkinter.CTkEntry(
            self.input_frame, placeholder_text='REF allele', font=('Helvetica', 30, 'bold')
        )
        self.ref.grid(row=row, column=1, pady=12, padx=10, sticky='ew')
        self.alt = customtkinter.CTkEntry(
            self.input_frame, placeholder_text='ALT allele', font=('Helvetica', 30, 'bold')
        )
        self.alt.grid(row=row + 1, column=0, pady=12, padx=10, sticky='ew')
        self.button = customtkinter.CTkButton(
            self.input_frame, text='Run', command=self.run_prime, font=('Helvetica', 50, 'bold')
        )
        self.button.grid(row=row + 2, column=0, pady=12, padx=10, sticky='ew', columnspan=2)
        self.output_frame = customtkinter.CTkFrame(self)
        self.output_frame.grid(row=2, column=0, padx=0, pady=(10, 0), sticky='ew')
        self.output1 = customtkinter.CTkTextbox(self.output_frame)
        self.output1.grid(row=0, column=0, pady=12, padx=10, sticky='ew', columnspan=3)
        self.button_3 = customtkinter.CTkButton(
            self.output_frame, text='Find new gRNAs', command=self.reinstate_gRNA_elements, font=('Helvetica', 50, 'bold')
        )
        self.button_3.grid(row=1, column=0, pady=12, padx=10, sticky='ew', columnspan=2)

    def run_prime(self):
        seq = self.entry1.get('0.0', 'end').strip()
        candidates = design_prime_edit(
            seq,
            PrimeEdit(int(self.position.get()), self.ref.get(), self.alt.get()),
            resolve_nuclease(),
        )
        lines = [
            f"#{c.rank} {c.spacer.spacer} PBS={c.pbs_length} RTT={c.rtt_length} "
            f"extension={c.extension_sequence} score={c.scores['structural_prime']}"
            for c in candidates[:200]
        ]
        self.output1.delete('0.0', 'end')
        self.output1.insert('0.0', '\n'.join(lines))


def main() -> None:
    app = App()
    app.mainloop()


if __name__ == '__main__':
    main()
