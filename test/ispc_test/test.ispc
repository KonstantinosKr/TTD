export void ispc_test (uniform int n, uniform int m,
                           uniform REAL * uniform a[3], uniform REAL * uniform b[3], uniform REAL * uniform c[3],
                           uniform REAL * uniform d[3], uniform REAL * uniform e[3], uniform REAL * uniform f[3],
                           uniform REAL * uniform p[3], uniform REAL * uniform q[3], uniform REAL distance[])
{
  foreach (i = 0 ... n)
  {
    REAL aa[3], bb[3], cc[3], dd[3], ee[3], ff[3];

    aa[0] = a[0][i];
    aa[1] = a[1][i];
    aa[2] = a[2][i];

    bb[0] = b[0][i];
    bb[1] = b[1][i];
    bb[2] = b[2][i];

    cc[0] = c[0][i];
    cc[1] = c[1][i];
    cc[2] = c[2][i];

    dd[0] = d[0][i];
    dd[1] = d[1][i];
    dd[2] = d[2][i];

    ee[0] = e[0][i];
    ee[1] = e[1][i];
    ee[2] = e[2][i];

    ff[0] = f[0][i];
    ff[1] = f[1][i];
    ff[2] = f[2][i];


    REAL out;
    for (int j = 0; j < m; j ++)
    {
      out = (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])+ 
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2]);
    }

    p[0][i] = out;
    p[1][i] = out;
    p[2][i] = out;

    q[0][i] = out;
    q[1][i] = out;
    q[2][i] = out;
  }
}

export uniform REAL * uniform ispc_allocate (uniform int size)
{
    return uniform new uniform REAL [size];
}
