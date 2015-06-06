table bonus_bed_columns "bonus columns"
(
    string chrom; "Chromosome"
    uint chromStart; "chr start"
    uint chromEnd; "chr end"
    string name; "item name"
    float my_floats; "some float values"
    set(item1,item2,item3,item4) my_sets; "some set options"
    int     my_ints; "signed integer values"
    lstring my_strs; "str representation of transcripts"
    uint[3] my_colors; "r,g,b colors"
)
