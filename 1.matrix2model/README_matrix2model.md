## Files (data, code and output) and directory tree:
 
````````````
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
1.matrix2model/ #extract draft models using "MEMs" (iMAT, GIMME, MBA, FastCore)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    1.Zhang/

		(#phenotype- primary astrocytes)
		(#v12 implements the method described in manuscript)
		
		* v12/abs/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_abs.m (MEMs_fpkm_abs) #code
			* matrix2models_abs_v12.mat (Models_fpkm_abs) #output
		* v12/norm_t1/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_norm_t1.m (MEMs_fpkm_norm_t1) #code
			* matrix2models_norm_t1_v12.mat (Models_fpkm_norm_t1) #output	
		* v12/norm_t2/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_norm_t2.m (MEMs_fpkm_norm_t2) #code
			* matrix2models_norm_t2_v12.mat (Models_fpkm_norm_t1) #output
			
    2.Vadodaria/

		(#phenotype- iPS-derived astrocytes from BD patients and controls)
		
		(#Below mentioned only for "1.Control_Untreated". In similar logic, 
		data & codes for "2.BD_Untreated", "3.BD_Responder_Untreated" 
		and "4.BD_NonResponder_Untreated" are available under "2.Vadodaria/")
				
		(#v3 implements the method described in manuscript)
		
		* 1.Control_Untreated/v3/abs/
			* Vadodaria_Control_Untreated.mat (ExpressionMatrix) #data
			* matrix2models_abs_vadodaria.m (MEMs_fpkm_abs) #code
			* matrix2models_abs_v3.mat (Models_fpkm_abs) #output
		* 1.Control_Untreated/v3/norm_t1/
			* Vadodaria_Control_Untreated.mat (ExpressionMatrix) #data
			* matrix2models_norm_t1_vadodaria.m (MEMs_fpkm_norm_t1) #code
			* matrix2models_norm_t1_v3.mat (Models_fpkm_norm_t1) #output
		* 1.Control_Untreated/v3/norm_t2/
			* Vadodaria_Control_Untreated.mat (ExpressionMatrix) #data
			* matrix2models_norm_t2_vadodaria.m (MEMs_fpkm_norm_t2) #code
			* matrix2models_norm_t2_v3.mat (Models_fpkm_norm_t2) #output
    
    3.Koskuvi/

		(#phenotype- iPS-derived astrocytes from monozygotic twin pairs discordant 
		for schizophrenia and healthy subjects)
		
		(#Below mentioned only for "1.Control"  (healthy controls). In similar logic, 
		data & codes for "2.HT" (healthy-twin) and "3.ST" (schizophrenia-twin) are 
		available under "3.Koskuvi/")

		(#v1 implements the method described in manuscript)

		* 1.Control/v1/abs/
			* Koskuvi_Control.mat (ExpressionMatrix) #data
			* matrix2models_abs_koskuvi.m (MEMs_fpkm_abs) #code
			* matrix2models_abs_v1.mat (Models_fpkm_abs) #output
		* 1.Control/v1/norm_t1/
			* Koskuvi_Control.mat (ExpressionMatrix) #data
			* matrix2models_norm_t1_koskuvi.m (MEMs_fpkm_norm_t1) #code
			* matrix2models_norm_t1_v1.mat (Models_fpkm_norm_t1) #output
		* 1.Control/v1/norm_t2/
			* Koskuvi_Control.mat (ExpressionMatrix) #data
			* matrix2models_norm_t2_koskuvi.m (MEMs_fpkm_norm_t2) #code
			* matrix2models_norm_t2_v1.mat (Models_fpkm_norm_t2) #output

````````````


