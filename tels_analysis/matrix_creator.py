import csv
import os

class matrix_creator:

    file_path = lambda self, file_name, extension : self.source_prefix + file_name + self.source_suffix + extension

    def __init__(self, source_prefix, source_suffix, amr_ext, mge_ext):
        self.source_prefix = source_prefix
        self.source_suffix = source_suffix
        self.amr_ext = amr_ext
        self.mge_ext = mge_ext
        self.amr_rows = dict()
        self.mge_rows = dict()
        self.column_names = list()

    def add_column(self, file_name):
        amr_file = self.file_path(file_name, self.amr_ext)
        mge_file = self.file_path(file_name, self.mge_ext)
        if os.path.exists(amr_file):
            with open(amr_file, 'r') as csv_file:
                csv_reader = csv.reader(csv_file)
                reader_row_num = 0
                for reader_row in csv_reader:
                    reader_row_num += 1
                    if reader_row_num < 20:
                        continue
                    if reader_row[0] not in self.amr_rows:
                        self.amr_rows[reader_row[0]] = ['0'] * len(self.column_names)
                    self.amr_rows[reader_row[0]].append(reader_row[1])

            with open(mge_file, 'r') as csv_file:
                csv_reader = csv.reader(csv_file)
                reader_row_num = 0
                for reader_row in csv_reader:
                    reader_row_num += 1
                    if reader_row_num < 20:
                        continue
                    if reader_row[0] not in self.mge_rows:
                        self.mge_rows[reader_row[0]] = ['0'] * len(self.column_names)
                    self.mge_rows[reader_row[0]].append(reader_row[1])

            self.column_names.append(file_name)
            for row in self.amr_rows:
                if len(self.amr_rows[row]) < len(self.column_names):
                    self.amr_rows[row].append('0')
            for row in self.mge_rows:
                if len(self.mge_rows[row]) < len(self.column_names):
                    self.mge_rows[row].append('0')

    def print_matrix(self, amr_matrix, mge_matrix):
        with open(amr_matrix, 'w') as csv_file:
            for name in self.column_names:
                csv_file.write(',' + name)
            csv_file.write('\n')
            for row_name, row_list in self.amr_rows.items():
                csv_file.write(row_name)
                for count in row_list:
                    csv_file.write(',' + count)
                csv_file.write('\n')

        with open(mge_matrix, 'w') as csv_file:
            for name in self.column_names:
                csv_file.write(',' + name)
            csv_file.write('\n')
            for row_name, row_list in self.mge_rows.items():
                csv_file.write(row_name)
                for count in row_list:
                    csv_file.write(',' + count)
                csv_file.write('\n')
