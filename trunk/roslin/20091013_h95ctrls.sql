create table "FANTOM4".h95ctrls (
  id varchar,
  cd14 int,
  control int,
  total int
  );
comment on table h95ctrls is 'Human cage tags from stimulated and control macroph. Source: fantom4 h95ctrls_tag_counts.tbl.txt.bz2';
comment on column h95ctrls.cd14 is 'Tag count in stimulated cells';
comment on column h95ctrls.control is 'Tag count in control cells';

-- Comment and header lines removed from h95ctrls_tag_counts.tbl.txt before executing COPY.
copy h95ctrls from E'C:\\Documents and Settings\\postgres\\My Documents\\Tritume\\h95ctrls_tag_counts.tbl.txt'
with delimiter E'\t'

alter table h95ctrls add constraint "PK_h95ctrls.id" primary key (id);