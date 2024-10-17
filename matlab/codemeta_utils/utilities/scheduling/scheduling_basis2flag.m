function flag=scheduling_basis2flag(basis,nbDays)
flag=reshape(basis,[],nbDays)'>0.5;
