
FPType genHH0(int i,ising1D* system)
{
  return system->get_h();
}


FPType genJJ0(int i,ising1D* system)
{
  return system->get_J();
}


FPType genHH(int i, ising1D* system)
{
  int size = system->get_size();

  if (( i >= size/2 - delta )&& (i < size/2 + delta)  )
    return system->get_h()+system->get_epsilon();
  else
    return system->get_h();
}


FPType genJJ(int i,ising1D* system)
{
  return system->get_J();
}



FPType genHHmod(int i, ising1D* system)
{
  int size = system->get_size();

  if ( ( (i/delta) % 2 ) ==0)
    return system->get_h()+system->get_epsilon();
  else
    return system->get_h();
}

FPType genHHpeaks(int i, ising1D* system)
{
  int size = system->get_size();

  if ( (i % delta) ==0)
    return system->get_h()+system->get_epsilon();
  else
    return system->get_h();
}
