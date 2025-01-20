drop table if exists af_structure cascade;
drop table if exists af_pdb_hit cascade;
drop table if exists af_transplant cascade;

create table af_structure (
	id serial primary key,
    name varchar not null,
	chunked boolean default false,
	af_version varchar not null,
    created timestamp with time zone not null,
	af_file varchar not null
);

alter table af_structure owner to "$OWNER";

create table af_pdb_hit (
	id serial primary key,
	af_id bigint references af_structure on delete cascade deferrable initially deferred,
	identity float,
	length integer,
	pdb_asym_id varchar not null,
	pdb_id varchar not null,
	rmsd float
);

alter table af_pdb_hit owner to "$OWNER";

create table af_transplant (
	id serial primary key,
	hit_id bigint references af_pdb_hit on delete cascade deferrable initially deferred,
	asym_id varchar not null,
	compound_id varchar not null,
	analogue_id varchar,
	entity_id varchar,
	rmsd float
);

alter table af_transplant owner to "$OWNER";

copy af_structure (id, name, chunked, af_version, created, af_file) from stdin;
$AF_STRUCTURES\.

copy af_pdb_hit (id, af_id, identity, length, pdb_asym_id, pdb_id, rmsd) from stdin;
$AF_PDB_HITS\.

copy af_transplant (id, hit_id, asym_id, compound_id, analogue_id, entity_id, rmsd) from stdin;
$AF_TRANSPLANTS\.

-- indices

create index hit_identity_ix on af_pdb_hit(identity);
create index hit_af_id_ix on af_pdb_hit(af_id);
