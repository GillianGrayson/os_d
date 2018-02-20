#include "processor.h"

Processor::Processor(AllData * ad)
{
	this->ad = ad;
}

void Processor::set_debug_behaviour(DebugBehavior* db)
{
	this->db = db;
}

void Processor::set_output_behaviour(OuputBehavior* ob)
{
	this->ob = ob;
}

void Processor::set_init_behaviour(NewDelBehavior* ndb)
{
	this->ndb = ndb;
}

void Processor::set_init_behaviour_qj(QJNewDelBehavior* ndb_qj)
{
	this->ndb_qj = ndb_qj;
}

void Processor::set_experiment_behaviour(ExperimentBehavior* eb)
{
	this->eb = eb;
}

void Processor::process()
{
	ndb->init_sizes(rp, cp, md);
	ndb->init_hamiltonians(rp, cp, md);
	ndb->init_dissipators(rp, cp, md);
	ndb->init_hamiltonians_qj(rp, cp, md);

	ob->suffix_param(cp, 4);
	db->save(rp, cp, md);

	ndb_qj->init_data(rp, cp, md, qjd);

	eb->trans_process(rp, cp, md, qjd);
	eb->obser_process(rp, cp, md, qjd);

	ndb_qj->free_data(rp, cp, md, qjd);

	ndb->free_hamiltonians(rp, cp, md);
	ndb->free_dissipators(rp, cp, md);
	ndb->free_hamiltonians_qj(rp, cp, md);
}