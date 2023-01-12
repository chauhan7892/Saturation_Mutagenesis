class MappedResult:
    # def __init__(self, seq_id, coord_position, seq_len, seq):
    def __init__(self, seq_id, coord_position, seq_len):
        self.seq_id = seq_id
        self.coord_position_first = coord_position
        self.seq_len = seq_len
        # self.seq = seq
