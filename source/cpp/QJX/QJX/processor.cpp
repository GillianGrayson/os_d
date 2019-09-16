#include "processor.h"

Processor::Processor(AllData * ad)
{
	this->ad = ad;
	this->db = nullptr;
	this->ob = nullptr;
	this->ndb = nullptr;
	this->ndb_exp = nullptr;
	this->eb = nullptr;
	this->pb = nullptr;
	this->cb = nullptr;
}

void Processor::set_debug_behaviour(DebugBehavior* db)
{
	this->db = db;
}

void Processor::set_output_behaviour(OutputBehavior* ob)
{
	this->ob = ob;
}

void Processor::set_newdel_behaviour(NewDelBehavior* ndb)
{
	this->ndb = ndb;
}

void Processor::set_newdel_exp_behaviour(ExpNewDelBehavior* ndb_exp)
{
	this->ndb_exp = ndb_exp;
}

void Processor::set_experiment_behaviour(ExperimentBehavior* eb)
{
	this->eb = eb;
}

void Processor::set_propagate_behaviour(PropagateBehavior* pb)
{
	this->pb = pb;
}

void Processor::set_core_behaviour(CoreBehavior* cb)
{
	this->cb = cb;
}

void Processor::process()
{
	ob->suffix_param(ad->rp, ad->cp, 4);

	ndb->init_sizes(ad);
	cout << "init_sizes" << endl;
	ndb->init_hamiltonians(ad);
	cout << "init_hamiltonians" << endl;
	ndb->init_dissipators(ad);
	cout << "init_dissipators" << endl;
	ndb->init_hamiltonians_qj(ad);
	cout << "init_hamiltonians_qj" << endl;

	db->save(ad);

	ndb_exp->init_data(ad, pb, cb);
	cout << "init_data" << endl;

	eb->trans_process(ad, pb, cb);
	eb->obser_process(ad, pb, cb);

	ndb_exp->free_data(ad, pb, cb);

	ndb->free_hamiltonians(ad);
	ndb->free_dissipators(ad);
	ndb->free_hamiltonians_qj(ad);
}