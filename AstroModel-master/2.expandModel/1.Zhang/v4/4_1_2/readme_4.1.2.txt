Get from 'MBA_model(v4.1)' to 'MBA_model_Synapse(v4.1.2)':

**********************************************************************

>> load('v4.1.mat')
>> open /media/anirudh/Work/ADBS_NIMHANS/Thesis/AstroModel/Databases_Tools/cobratoolbox/AstroModel_QC/2. sanityCheck/4.1.1/v4.1.1_inter2.xlsx
>> create addSynapseRxns.m based on "v4.1.1_inter2.xlsx" & "addReaction()"
>> save & run addSynapseRxns.m

>> copy metabolite information (charge, formula) from "v4.1.1_inter2.xlsx" to the ("metCharges" & "metFormulas") fields in MBA_model
>> save MBA_model_Synapse as "v4.1.2.mat"

**********************************************************************

>> load model v4.1.2.mat
>> add glucose exchange {'EX_glc_D[e]'} to model
>> remove glucose_transport synapse reaction {'GLCt4_SYN'} from model
>> model size: [5403 x 9627]
>> remove oxygen exchange synapse reaction {'EX_o2[s]'} from model
>> model size: [5402 x 9626]
>> add oxygen transport diffusion {'O2t'} to model
>> model size: [5402 x 9627]
>> add ATP demand reaction {'DM_atp_c_'} to model
>> model size: [5402 x 9628]

>> Set 'EX_o2[e]' 'ub' to zero. Normally cells don’t excrete oxygen. It’s only uptaken

**********************************************************************
>> #23 Rxns involve 'k[c]' in present astrocyte metabolic model.. 
>> #12 Rxns involve 'k[e]' in present astrocyte metabolic model.. 
**********************************************************************

Test4HumanBrain.m:
BrainTasks_old:
>> add ATP demand reaction {'EX_adn[e]','sink_amp[c]','sink_imp[c]','sink_gmp[c]','sink_melatn[c]','sink_hista[c]'} to model
>> model size: [5402 x 9634]
BrainTasks_New:



**********************************************************************
Model stats (05-Nov):
Model size: [5402 x 9634]
Unique mets: 2312
Unique genes: 2052
Compartments: 10
*********************
Rxns	Blood	Synapse
Exchange (N=1607)	805	802	%blood exclusive {'EX_adn', 'EX_o2', 'EX_glc_D'}%
Demand (N=50)	-	-
Sink (N=51)	-	-
**********************************************************************
>> model size: [5402 x 9634]
>> objective: 'biomass_maintenance'
>> FBA: model.f = 215.36
>> only 3357 Rxns carry a ~0 flux, for this objective.
**********************************************************************


