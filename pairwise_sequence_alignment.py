import sys, os

class ScoringMatrix:
    def __init__(self, matrix, row_titles, column_titles, row_index_map, column_index_map):
        self.matrix = matrix
        self.row_titles = row_titles
        self.column_titles = column_titles
        self.row_index_map = row_index_map
        self.column_index_map = column_index_map
    
    def get_matrix(self):
        return self.matrix
    
    def get_row_titles(self):
        return self.row_titles
    
    def get_column_titles(self):
        return self.column_titles

    def get_row_index_map(self):
        return self.row_index_map
    
    def get_column_index_map(self):
        return self.column_index_map
    
    def get_score_of_pair(self, x, y):
        row_index = self.row_index_map[x.upper()]
        column_index = self.column_index_map[y.upper()]
        return self.matrix[row_index][column_index]
    
    def print(self):
        max_len = len(str(self.matrix[0][0]))
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                l = len(str(self.matrix[i][j]))
                if l > max_len:
                    max_len = l
        
        space = " " * max_len

        print (" " + space, end = '')
        print ((space).join(self.column_titles))
        for i in range(len(self.matrix)):
            print (self.row_titles[i], end = ' ')
            for j in range(len(self.matrix[0])):
                print (str(self.matrix[i][j]).rjust(max_len), end = ' ')
            print ()

class ScoringMatrixFileReader:
    def __init__(self, path = None):
        self.path = path

    def load_path(self, path):
        self.path = path

    def read_matrix(self):
        valid_lines = []
        with open(self.path, encoding="utf-8") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if not line.strip().startswith("#"):
                    valid_lines.append(line.strip().split("#")[0].upper())

        column_titles = valid_lines[0].split()
        column_index_map = {i:column_titles.index(i) for i in column_titles}

        row_titles = [i.split()[0] for i in valid_lines[1:]]
        row_index_map = {i:row_titles.index(i) for i in row_titles}

        matrix = []
        for line in valid_lines[1:]:
            values = []
            for value in line.split()[1:]:
                values.append(int(value))
            matrix.append(values)

        return ScoringMatrix(matrix, row_titles, column_titles, row_index_map, column_index_map)

class AlignmentProcessor:
    # the algorithm

    def __init__(self, sequence1, sequence2, alignment_type,
     scoring_matrix, gap_opening_penalty, gap_extension_penalty, DEBUG):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.sequence1_length = len(self.sequence1)
        self.sequence2_length = len(self.sequence2)
        self.alignment_type = alignment_type
        self.scoring_matrix = scoring_matrix
        self.gap_opening_penalty = gap_opening_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.algorithm_matrix = []
        self.DEBUG = DEBUG
    
    def calculate_percent_identities(self, list_of_sequence_pairs, round_to = 4):
        results = []
        for sequence1, sequence2 in list_of_sequence_pairs:
            number_of_matches = 0
            n = len(sequence1)
            for i in range(n):
                if sequence1[i] == sequence2[i]:
                    number_of_matches += 1
            
            result = number_of_matches * 100 / n
            factor = 10 ** round_to
            results.append(round(result * factor) / factor)

        return results


    def calculate_raw_alignment_scores(self, list_of_sequence_pairs):
        results = []

        for sequence1, sequence2 in list_of_sequence_pairs:
            nongap_score = 0
            opening_gaps = 0
            extended_gaps = 0

            gap_extension_state_seq1 = False
            gap_extension_state_seq2 = False

            for i in range(len(sequence1)):
                from_seq1 = sequence1[i]
                from_seq2 = sequence2[i]

                if from_seq1 != '-' and from_seq2 != '-':
                    if gap_extension_state_seq1:
                        gap_extension_state_seq1 = False
                    if gap_extension_state_seq2:
                        gap_extension_state_seq2 = False

                    nongap_score += self.scoring_matrix.get_score_of_pair(from_seq1, from_seq2)

                else:
                    if from_seq1 == '-':
                        if gap_extension_state_seq1:
                            extended_gaps += 1
                        else:
                            opening_gaps += 1
                            gap_extension_state_seq1 = True
                    
                    if from_seq2 == '-':
                        if gap_extension_state_seq2:
                            extended_gaps += 1
                        else:
                            opening_gaps += 1
                            gap_extension_state_seq2 = True
        
            results.append(nongap_score + (opening_gaps * self.gap_opening_penalty) + (extended_gaps * self.gap_extension_penalty))
        return results


    def align(self):
        if self.alignment_type == 'global':
            algorithm = self.global_alignment
        elif self.alignment_type == 'local':
            algorithm = self.local_alignment

        results = algorithm()
        sequence_match_strings = self.get_sequence_report_strings(results)
        raw_alignment_scores = self.calculate_raw_alignment_scores(results)
        percent_identities = self.calculate_percent_identities(results)
        return results, sequence_match_strings, raw_alignment_scores, percent_identities

    def create_algorithm_matrix(self):
        if self.alignment_type == 'global':
            m, n = self.create_algorithm_matrix_for_global_alignment()
        else:
            m, n = self.create_algorithm_matrix_for_local_alignment()
        
        return m, n
    
    def create_algorith_matrix_skeleton(self):
        m = self.sequence2_length + 1 # number of rows
        n = self.sequence1_length + 1 # number of columns

        self.algorithm_matrix = [[0 for j in range(n)] for i in range(m)] # m x n matrix

        return m, n

    def create_algorithm_matrix_for_local_alignment(self):
        m, n = self.create_algorith_matrix_skeleton()

        for i in range(m):
            self.algorithm_matrix[i][0] = 0
        for j in range(n):
            self.algorithm_matrix[0][j] = 0
        
        return m, n

    def create_algorithm_matrix_for_global_alignment(self):
        m, n = self.create_algorith_matrix_skeleton()
        self.algorithm_matrix[0][0] = 0

        for i in range(1, min(2, m)):
            self.algorithm_matrix[i][0] = self.gap_opening_penalty
        for j in range(1, min(2, n)):
            self.algorithm_matrix[0][j] = self.gap_opening_penalty

        for i in range(2, m):
            self.algorithm_matrix[i][0] = self.algorithm_matrix[i - 1][0] + self.gap_extension_penalty
        for j in range(2, n):
            self.algorithm_matrix[0][j] = self.algorithm_matrix[0][j - 1] + self.gap_extension_penalty
        
        return m, n

    def get_sequence_report_strings(self, list_of_sequence_pairs):
        results = []
        for x in list_of_sequence_pairs:
            s = ""
            seq1 = x[0]
            seq2 = x[1]
            for i in range(len(seq1)):
                if seq1[i] == seq2[i]:
                    s += '|'
                else:
                    s += ' '
            results.append(s)
        return results

    def local_alignment(self):
        # -1: gap in sequence1 (vertical movement)
        #  0: diagonal movement
        #  1: gap in sequence2 (horizontal movement)
        #  2: 0 (none of above)

        movement_map = {}
        m, n = self.create_algorithm_matrix()

        for i in range(1, m):
            for j in range(1, n):
                aminoacid_from_sequence1 = self.sequence1[j - 1]
                aminoacid_from_sequence2 = self.sequence2[i - 1]

                if (i - 1, j) in movement_map and movement_map[(i - 1, j)] == -1:
                    vertical_score = self.gap_extension_penalty
                else:
                    vertical_score = self.gap_opening_penalty
                
                if (i, j - 1) in movement_map and movement_map[(i, j - 1)] == 1:
                    horizontal_score = self.gap_extension_penalty
                else:
                    horizontal_score = self.gap_opening_penalty
                
                diagonal_score = self.scoring_matrix.get_score_of_pair(aminoacid_from_sequence1, aminoacid_from_sequence2)

                vertical_movement_consequence = self.algorithm_matrix[i - 1][j] + vertical_score
                diagonal_movement_consequence = self.algorithm_matrix[i - 1][j - 1] + diagonal_score
                horizontal_movement_consequence = self.algorithm_matrix[i][j - 1] + horizontal_score

                maximum = max(vertical_movement_consequence, diagonal_movement_consequence, horizontal_movement_consequence)

                if maximum < 0:
                    self.algorithm_matrix[i][j] = 0
                    movement_map[(i, j)] = 2
                else:
                    self.algorithm_matrix[i][j] = maximum
                    if maximum == diagonal_movement_consequence:
                        movement_map[(i, j)] = 0
                    elif maximum == vertical_movement_consequence:
                        movement_map[(i, j)] = -1
                    else:
                        movement_map[(i, j)] = 1
                    
        if self.DEBUG:
            for i in range(m):
                for j in range(n):
                    print (str(self.algorithm_matrix[i][j]).rjust(4), end= ' ')
                print ()

        # go from the bottom to the top: compose the alignment
        top_line = "" # aligned sequence1
        bottom_line = "" # aligned sequence2

        max_locations = self.find_max_locations_in_the_matrix(self.algorithm_matrix)
        if self.DEBUG:
            print (max_locations)


        results = []
        for max_location in max_locations:
            i, j = max_location

            while self.algorithm_matrix[i][j] > 0 and i > 0 and j > 0:
                movement = movement_map[(i, j)]
                if movement == 0:
                    top_line = self.sequence1[j - 1] + top_line
                    bottom_line = self.sequence2[i - 1] + bottom_line
                    i -= 1
                    j -= 1
                elif movement == -1:
                    top_line = "-" + top_line
                    bottom_line = self.sequence2[i - 1] + bottom_line
                    i -= 1
                elif movement == 1:
                    top_line = self.sequence1[j - 1] + top_line
                    bottom_line = "-" + bottom_line
                    j -= 1

            if self.DEBUG:
                print(movement_map)

            results.append((top_line, bottom_line))

        return results
    
    def find_max_locations_in_the_matrix(self, matrix):
        # matrix is a m x n matrix.
        m = len(matrix)
        n = len(matrix[0])

        max_val = matrix[0][0]

        for i in range(m):
            for j in range(n):
                if matrix[i][j] > max_val:
                    max_val = matrix[i][j]

        max_locations = []

        for i in range(m):
            for j in range(n):
                if matrix[i][j] == max_val:
                    max_locations.append((i, j))

        return max_locations

    def global_alignment(self):
        # -1: gap in sequence1 (vertical movement)
        #  0: diagonal movement
        #  1: gap in sequence2 (horizontal movement)

        movement_map = {}

        m, n = self.create_algorithm_matrix()
        # fill algorithm matrix
        for i in range(1, m):
            for j in range(1, n):
                aminoacid_from_sequence1 = self.sequence1[j - 1]
                aminoacid_from_sequence2 = self.sequence2[i - 1]
                if (i - 1, j) in movement_map and movement_map[(i - 1, j)] == -1:
                    vertical_score = self.gap_extension_penalty
                else:
                    vertical_score = self.gap_opening_penalty

                diagonal_score = self.scoring_matrix.get_score_of_pair(aminoacid_from_sequence1, aminoacid_from_sequence2)

                if (i, j - 1) in movement_map and movement_map[(i, j - 1)] == 1:
                    horizontal_score = self.gap_extension_penalty
                else:
                    horizontal_score = self.gap_opening_penalty

                vertical_movement_consequence = self.algorithm_matrix[i - 1][j] + vertical_score
                diagonal_movement_consequence = self.algorithm_matrix[i - 1][j - 1] + diagonal_score
                horizontal_movement_consequence = self.algorithm_matrix[i][j - 1] + horizontal_score

                maximum = max(vertical_movement_consequence, diagonal_movement_consequence, horizontal_movement_consequence)
                self.algorithm_matrix[i][j] = maximum
                if maximum == diagonal_movement_consequence:
                    movement_map[(i, j)] = 0
                elif maximum == vertical_movement_consequence:
                    movement_map[(i, j)] = -1
                else:
                    movement_map[(i, j)] = 1
        
        if self.DEBUG:
            print (movement_map)

        # go from the bottom to the top: compose the alignment
        top_line = "" # aligned sequence1
        bottom_line = "" # aligned sequence2

        i = m - 1
        j = n - 1
        while i > 0 and j > 0:
            movement = movement_map[(i, j)]
            if movement == 0:
                top_line = self.sequence1[j - 1] + top_line
                bottom_line = self.sequence2[i - 1] + bottom_line
                i -= 1
                j -= 1
            elif movement == -1:
                top_line = "-" + top_line
                bottom_line = self.sequence2[i - 1] + bottom_line
                i -= 1
            else:
                top_line = self.sequence1[j - 1] + top_line
                bottom_line = "-" + bottom_line
                j -= 1

        while i > 0:
            top_line = "-" + top_line
            bottom_line = self.sequence2[i - 1] + bottom_line
            i -= 1
        
        while j > 0:
            top_line = self.sequence1[j - 1] + top_line
            bottom_line = "-" + bottom_line
            j -= 1

        if self.DEBUG:
            for i in range(m):
                for j in range(n):
                    print (str(self.algorithm_matrix[i][j]).rjust(4), end= ' ')
                print ()

        return [(top_line, bottom_line)]

class Main:
    def __init__(self, args, DEBUG):
        self.DEBUG = DEBUG
        self.args = args
        input_path, alignment_type, scoring_matrix_path, gap_opening_penalty, gap_extension_penalty, output_path = self.check_args(args)
        initial_error = False
        output_file = False

        numbers = [gap_opening_penalty, gap_extension_penalty]
        number_names = ('gap opening penalty', 'gap extension penalty')
        
        for i in range(len(numbers)):
            temp_error = False
            name = number_names[i]
            try:
                numbers[i] = int(numbers[i])
            except ValueError:
                print ("\nInvalid value for {}!: '{}'".format(name, numbers[i]))
                initial_error = True
                temp_error = True
            number = numbers[i]
            if not temp_error and number > 0:
                print ("\n{} cannot be positive!".format(name.capitalize()))
                initial_error = True
        gap_opening_penalty, gap_extension_penalty = numbers

        paths = [input_path, scoring_matrix_path]
        path_names = ('input path', 'scoring matrix path')
        for i in range(len(paths)):
            name = path_names[i]
            paths[i] = os.path.abspath(paths[i])
            path = paths[i]
            if not os.path.exists(path):
                print ("\nThe path given for {} does not exist!: {}".format(name, path))
                initial_error = True
        input_path, scoring_matrix_path = paths

        if alignment_type not in ('local', 'global'):
            print ('\nInvalid value for alignment type!: {}\nAlignment type can be either local or global'.format(alignment_type))
            initial_error = True

        if output_path is not None and os.path.exists(output_path):
            print ("\nThe specified output path already exists!: '{}'\nPlease specify a non-existing output path.".format(output_path))
            initial_error = True
        
        if initial_error:
            sys.exit()

        if output_path is not None:
            try:        
                with open(output_path, 'w') as f:
                    pass
            except Exception as e:
                print ("Cannot use the path {} for writing!\nError!: {}\nPlease specify a valid path for writing.".format(output_path, str(e)))
            output_file = True 
            output_path = os.path.abspath(output_path)           


        sequence1, sequence2 = self.read_input(input_path)
        if self.DEBUG:
            print (sequence1, sequence2)

        score_matrix_file_reader = ScoringMatrixFileReader(scoring_matrix_path)
        scoring_matrix = score_matrix_file_reader.read_matrix()

        alignment_processor = AlignmentProcessor(sequence1, sequence2, alignment_type,
                                                    scoring_matrix, gap_opening_penalty,
                                                    gap_extension_penalty, self.DEBUG)

        aligned_sequences, match_strings, raw_alignment_scores, percent_identities = alignment_processor.align()

        # The first part of the output
        if not output_file:
            for i in range(len(aligned_sequences)):
                print (aligned_sequences[i][0])
                print (match_strings[i])
                print (aligned_sequences[i][1])
                print ("Raw alignment score: {}".format(raw_alignment_scores[i]))
                print ("The percent identity between two aligned sequences: {}%".format(percent_identities[i]))

                if i != len(aligned_sequences) - 1:
                    print ("\n\n")

        else:
            with open(output_path, "w", encoding="utf-8") as f:
                for i in range(len(aligned_sequences)):
                    f.write(aligned_sequences[i][0] + "\n")
                    f.write(match_strings[i] + "\n")
                    f.write(aligned_sequences[i][1] + "\n")
                    if i != len(aligned_sequences) - 1:
                        f.write("\n\n")
            print ("The alignment output has been recorded to the following path: {}".format(output_path))
            if len(aligned_sequences) > 1:
                print ("There were more than 1 results. The results given below belong to the alignments recorded in the output file in the same order.")
            
            for i in range(len(aligned_sequences)):
                # The second part of the output
                print ("Raw alignment score: {}".format(raw_alignment_scores[i]))
                # The third part of the output
                print ("The percent identity between two aligned sequences: {}%".format(percent_identities[i]))

                if i != len(aligned_sequences) - 1:
                    print ("\n\n")

    def read_input(self, path):
        with open(path) as f:
            content = f.read().strip().splitlines()
        
        clear_content = [i for i in content if i.strip()]

        sequence1 = clear_content[0]
        sequence2 = clear_content[1]

        return (sequence1, sequence2)

    def check_args(self, args):
        if len(args) != 11 and len(args) != 13:
            self.print_usage_and_exit(args)

        expected_arg_markers = ("--input", "--alignment", "--scoring-matrix", "--gap-opening-penalty", "--gap-extension-penalty")
        optional_output_marker = "--output"
        arg_markers_indexes = []
        terminate = False
        for marker in expected_arg_markers:
            try:
                index = args.index(marker)
                arg_markers_indexes.append(index)
            except ValueError:
                terminate = True
                break

            if index > len(args) - 2 or args.count(marker) > 1:
                terminate = True
                break
        
        if terminate:
            self.print_usage_and_exit(args)

        if optional_output_marker in args:
            arg_markers_indexes.append(args.index(optional_output_marker))
            if arg_markers_indexes[-1] > len(args) - 2:
                self.print_usage_and_exit(args)
        else:
            arg_markers_indexes.append(None)

        results = []

        for i in range(len(arg_markers_indexes)):
            index = arg_markers_indexes[i]
            if index is None:
                results.append(None)
            else:
                results.append(args[index + 1])

        return results

    def print_usage(self, args):
        fn = os.path.split(args[0])[1]
        print ("Usage: python3 {} --input <path to input text file containing amino acid sequences> --alignment <local or global> --scoring-matrix <path to scoring matrix> --gap-opening-penalty <a negative number> --gap-extension-penalty <a negative number> --output <path to output file>\n--output is optional".format(fn))
        
    def print_usage_and_exit(self , args):
        self.print_usage(args)
        sys.exit()
        
        
if __name__ == "__main__":
    main = Main(sys.argv, DEBUG=False)