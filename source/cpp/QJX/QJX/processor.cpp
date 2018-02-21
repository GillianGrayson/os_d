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

void Processor::set_init_behaviour_qj(ExpNewDelBehavior* ndb_qj)
{
	this->ndb_qj = ndb_qj;
}

void Processor::set_experiment_behaviour(ExperimentBehavior* eb)
{
	this->eb = eb;
}

void Processor::set_propagate_behaviour(PropagateBehavior* pb)
{
	this->pb = pb;
}

void Processor::process()
{
	ndb->init_sizes(ad);
	ndb->init_hamiltonians(ad);
	ndb->init_dissipators(ad);
	ndb->init_hamiltonians_qj(ad);

	ob->suffix_param(ad->cp, 4);
	db->save(ad);

	ndb_qj->init_data(ad);

	eb->trans_process(ad, pb);
	eb->obser_process(ad, pb);

	ndb_qj->free_data(ad);

	ndb->free_hamiltonians(ad);
	ndb->free_dissipators(ad);
	ndb->free_hamiltonians_qj(ad);
}