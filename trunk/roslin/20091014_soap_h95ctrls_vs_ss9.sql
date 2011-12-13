CREATE TABLE "FANTOM4"."tmp_hg95ctrls_ss9"
(
  read_name character varying NOT NULL,
  read_sequence character varying(36) NOT NULL,
  read_quality character varying(36),
  "NoMatches" integer,
  "PairEnd" character varying(5),
  length_read integer,
  strand character varying(1),
  chromosome character varying(80),
  "position" integer,
  no_mismatches integer,
  "Ref_vs_Read" character varying(25)
);

copy "tmp_hg95ctrls_ss9" from E'C:\\Documents and Settings\\postgres\\My Documents\\Tritume\\20091014_hg95ctrls_vs_ss9_r1.mapr' with csv delimiter E'\t';