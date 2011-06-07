/*
Produce summary statistics for an alignment in SAM format produced by BWA
*/

select read_sam($$ file:'D:/Tritume/20101213_cage_vs_ss9_hg_rDNA_21bp.sam', dest_table:'tmp_sam_cage050810_21bp_bwa'  $$);
select read_sam($$ file:'D:/Tritume/20101213b_cage050810_iter22bp.sam', dest_table:'tmp_sam_cage050810_20101213b_iter22bp_bwa' $$);
select read_sam($$ file:'D:/Tritume/20101213b_cage050810_iter23bp.sam', dest_table:'tmp_sam_cage050810_20101213b_iter23bp_bwa' $$);
select read_sam($$ file:'D:/Tritume/20101213b_cage050810_iter24bp.sam', dest_table:'tmp_sam_cage050810_20101213b_iter24bp_bwa' $$);
select read_sam($$ file:'D:/Tritume/20101213b_cage050810_iter25bp.sam', dest_table:'tmp_sam_cage050810_20101213b_iter25bp_bwa' $$);

-- Unmapped
select * from (
    select 'Unmapped' AS category, count(qname) AS count from tmp_sam_cage050810_20101213b_iter25bp_bwa where flag & 4 = 4
    -- Mapped
    union select 'Mapped', count(qname) from tmp_sam_cage050810_20101213b_iter25bp_bwa where flag & 4 = 0
    -- hg rDNA
    union select 'Human rDNA', count(qname) from tmp_sam_cage050810_20101213b_iter25bp_bwa where rname like 'gi|555853|gb|U13369.1|HSU13369'
    -- Single and multi mappers
    union select sam_get_tagvalue(tags, 'X0:i'), count(qname) from tmp_sam_cage050810_20101213b_iter25bp_bwa 
      where sam_get_tagvalue(tags, 'X0:i')::int <= 10 
      group by sam_get_tagvalue(tags, 'X0:i') 
    -- Multimappers >10
    union select '10+', count(qname) from tmp_sam_cage050810_20101213b_iter25bp_bwa where sam_get_tagvalue(tags, 'X0:i')::int > 10
    ) as t1 order by category;
    