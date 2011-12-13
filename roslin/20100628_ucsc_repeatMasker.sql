/*
  Download and import repeat elements in the pig genome.  
  
  File 'genome_browser_pig_repeat_masker.txt' downloaded from USCS Genome Browser
  http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=163528763&clade=mammal&org=Pig&db=susScr2&hgta_group=varRep&hgta_track=susScr2&hgta_table=0&hgta_regionType=genome&position=chr13%3A57394166-57402412&hgta_outputType=primaryTable&hgta_outFileName=

*/

select read_table($$ file:'D:/Tritume/genome_browser_pig_repeat_masker.txt', dest_table: 'ucsc_repeat_masker_pig', header: True, limit:-1, overwrite:False $$)
select names($$ t: 'ucsc_repeat_masker_pig', rename: ['bin', 'sw_score', 'milli_div', 'milli_del', 'milli_ins', 'geno_name', 'geno_start', 'geno_end', 'geno_left', 'strand', 'rep_name', 'rep_class', 'rep_family', 'rep_start', 'rep_end', 'rep_left', 'id']$$)

COMMENT ON TABLE ucsc_repeat_masker_pig IS 'Table of repeat elements in the pig genome. Downloaded from UCSC on 28/06/2010. Pig genome version: SGSC Sscrofa9.2/susScr2'
COMMENT ON COLUMN ucsc_repeat_masker_pig.bin IS 'Indexing field to speed chromosome range queries.';
COMMENT ON COLUMN ucsc_repeat_masker_pig.sw_score IS 'Smith Waterman alignment score';
COMMENT ON COLUMN ucsc_repeat_masker_pig.milli_div IS 'Base mismatches in parts per thousand';
COMMENT ON COLUMN ucsc_repeat_masker_pig.milli_del IS 'Bases deleted in parts per thousand';
COMMENT ON COLUMN ucsc_repeat_masker_pig.milli_ins IS 'Bases inserted in parts per thousand';
COMMENT ON COLUMN ucsc_repeat_masker_pig.geno_name IS 'Genomic sequence name';
COMMENT ON COLUMN ucsc_repeat_masker_pig.geno_start IS 'Start in genomic sequence';
COMMENT ON COLUMN ucsc_repeat_masker_pig.geno_end IS 'End in genomic sequence';
COMMENT ON COLUMN ucsc_repeat_masker_pig.geno_left IS '-#bases after match in genomic sequence';
COMMENT ON COLUMN ucsc_repeat_masker_pig.strand IS 'Relative orientation + or -';
COMMENT ON COLUMN ucsc_repeat_masker_pig.rep_name IS 'Name of repeat';
COMMENT ON COLUMN ucsc_repeat_masker_pig.rep_class IS 'Class of repeat';
COMMENT ON COLUMN ucsc_repeat_masker_pig.rep_family IS 'Family of repeat';
COMMENT ON COLUMN ucsc_repeat_masker_pig.rep_start IS 'Start (if strand is +) or -#bases after match (if strand is -) in repeat sequence';
COMMENT ON COLUMN ucsc_repeat_masker_pig.rep_end IS 'End in repeat sequence';
COMMENT ON COLUMN ucsc_repeat_masker_pig.rep_left IS '-#bases after match (if strand is +) or start (if strand is -) in repeat sequence';
COMMENT ON COLUMN ucsc_repeat_masker_pig.id IS 'First digit of id field in RepeatMasker .out file. Best ignored.';

select * from ucsc_repeat_masker_pig limit 10;
select count(*) from ucsc_repeat_masker_pig ;  --4309043 rows

select rep_name, count(*) AS no_repeats from ucsc_repeat_masker_pig group by rep_name order by no_repeats;
 

