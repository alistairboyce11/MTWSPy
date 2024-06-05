class Phases:
    """
    Class to handle phases of interest for time delay picking
    
    """
    def __init__(self):
        pass

    def get_phase_dictionary(self):
        """
        :return output_dict: phases for each component to pick
        :type output_dict: dict
        """

        # Master phases
        CV_PHS_T = 's,S,SS,SSS,ScS,Sdiff,ScSScS,ScSScSScS,sS,sSS,sSSS,sScS,sSdiff,sScSScS,sScSScSScS'.split(',')
        CV_PHS_Z = 'p,P,PP,PPP,PcP,Pdiff,PKP,PKIKP,PcPPKP,PKIKPPKIKP,PcPPKPPKP,pP,sP,pPP,sPP,pPPP,sPPP,pPcP,sPcP,pPdiff,sPdiff,pPKP,sPKP,pPKIKP,sPKIKP,pPcPPKP,sPcPPKP,pPKIKPPKIKP,sPKIKPPKIKP,pPcPPKPPKP,sPcPPKPPKP,s,S,SS,SSS,ScS,Sdiff,SKS,SKIKS,ScP,SP,SPP,SSP,SKP,SKIKP,pS,sS,pSS,sSS,pSSS,sSSS,pScS,sScS,pSdiff,sSdiff,pSKS,sSKS,pScP,sScP,pSP,sSP,pSPP,sSPP,pSSP,sSSP,pSKP,sSKP'.split(',')
        CV_PHS_R = CV_PHS_Z

        # direct arrivals
        CV_PHS_T_1 = 's,S,Sdiff,ScS,sS,sSdiff,sScS'.split(',')
        CV_PHS_Z_1 = 'p,P,Pdiff,PcP,PKP,PKIKP,s,S,Sdiff,ScP,ScS,SKP,SKS,SKIKP,SKIKS,pP,sP,pPdiff,sPdiff,pPcP,sPcP,pPKP,sPKP,pPKIKP,sPKIKP,pS,sS,pSdiff,sSdiff,pScP,sScP,pSKP,sSKP,sSKS,pSKS,pSKIKP,sSKIKP,pSKIKS,sSKIKS'.split(',')
        CV_PHS_R_1 = CV_PHS_Z_1

        # surface reflections [1 bounce]
        CV_PHS_T_2 = CV_PHS_T_1 + 'SS,ScSScS,sSS,sScSScS'.split(',')
        CV_PHS_Z_2 = CV_PHS_Z_1 + 'PP,SP,SS,PcPPKP,PKPPKP,SKPPKP,SKPSKP,PKIKPPKIKP,pPP,sPP,pSP,sSP,pSS,sSS,pPcPPKP,sPcPPKP,pPKPPKP,sPKPPKP,pSKPPKP,sSKPPKP,pSKPSKP,sSKPSKP,pPKIKPPKIKP,sPKIKPPKIKP'.split(',')
        CV_PHS_R_2 = CV_PHS_Z_2

        # surface reflections [2 bounces]
        CV_PHS_T_3 = CV_PHS_T_2 + 'SSS,sSSS,ScSScSScS,sScSScSScS'.split(',')
        CV_PHS_Z_3 = CV_PHS_Z_2 + 'PPP,SPP,SSP,SSS,PcPPKPPKP,pPPP,sPPP,pSPP,sSPP,pSSP,sSSP,pSSS,sSSS,pPcPPKPPKP,sPcPPKPPKP'.split(',')
        CV_PHS_R_3 = CV_PHS_Z_3

        # surface reflections [3 bounces]
        CV_PHS_T_4 = CV_PHS_T_3 + 'SSSS,sSSSS'.split(',')
        CV_PHS_Z_4 = CV_PHS_Z_3 + 'PPPP,pPPPP,sPPPP,SSSS,pSSSS,sSSSS'.split(',')
        CV_PHS_R_4 = CV_PHS_Z_4

        # pure P
        CV_PHS_P_1 = 'p,P,Pdiff,PcP,PKP,PKIKP,pP,sP,pPdiff,pPcP,pPKP,pPKIKP'.split(',')
        CV_PHS_P_2 = CV_PHS_P_1 + 'PP,PcPPKP,PKPPKP,PKIKPPKIKP,pPP,pPcPPKP,pPKPPKP,pPKIKPPKIKP'.split(',')
        CV_PHS_P_3 = CV_PHS_P_2 + 'PPP,PcPPKPPKP,PKPPKPPKP,pPPP,pPcPPKPPKP,pPKPPKPPKP'.split(',')
        CV_PHS_P_4 = CV_PHS_P_3 + 'PPPP,pPPPP'.split(',')

        CV_PHASES = {'T': CV_PHS_T, 'Z': CV_PHS_Z, 'R': CV_PHS_R}
        CV_PHASES_1 = {'T': CV_PHS_T_1, 'Z': CV_PHS_Z_1, 'R': CV_PHS_R_1}
        CV_PHASES_2 = {'T': CV_PHS_T_2, 'Z': CV_PHS_Z_2, 'R': CV_PHS_R_2}
        CV_PHASES_3 = {'T': CV_PHS_T_3, 'Z': CV_PHS_Z_3, 'R': CV_PHS_R_3}
        CV_PHASES_4 = {'T': CV_PHS_T_4, 'Z': CV_PHS_Z_4, 'R': CV_PHS_R_4}

        output_dict = {'Master': CV_PHASES, 'Direct': CV_PHASES_1, '1_refl': CV_PHASES_2, '2_refl': CV_PHASES_3, '3_refl': CV_PHASES_4 }

        return output_dict


def main():
    phases_instance = Phases()
    output_dict = phases_instance.get_phase_dictionary()
    print(output_dict)


if __name__ == '__main__':
    main()







# # Common variables
# # Phasenames
# # For CV_PHASES(.Z|R|T) struct

# def phases():

#     # Master phases
#     CV_PHS_T = 's,S,SS,SSS,ScS,Sdiff,ScSScS,ScSScSScS,sS,sSS,sSSS,sScS,sSdiff,sScSScS,sScSScSScS'.split(',')
#     CV_PHS_Z = 'p,P,PP,PPP,PcP,Pdiff,PKP,PKIKP,PcPPKP,PKIKPPKIKP,PcPPKPPKP,pP,sP,pPP,sPP,pPPP,sPPP,pPcP,sPcP,pPdiff,sPdiff,pPKP,sPKP,pPKIKP,sPKIKP,pPcPPKP,sPcPPKP,pPKIKPPKIKP,sPKIKPPKIKP,pPcPPKPPKP,sPcPPKPPKP,s,S,SS,SSS,ScS,Sdiff,SKS,SKIKS,ScP,SP,SPP,SSP,SKP,SKIKP,pS,sS,pSS,sSS,pSSS,sSSS,pScS,sScS,pSdiff,sSdiff,pSKS,sSKS,pScP,sScP,pSP,sSP,pSPP,sSPP,pSSP,sSSP,pSKP,sSKP'.split(',')
#     CV_PHS_R = CV_PHS_Z

#     # direct arrivals
#     CV_PHS_T_1 = 's,S,Sdiff,ScS,sS,sSdiff,sScS'.split(',')
#     CV_PHS_Z_1 = 'p,P,Pdiff,PcP,PKP,PKIKP,s,S,Sdiff,ScP,ScS,SKP,SKS,SKIKP,SKIKS,pP,sP,pPdiff,sPdiff,pPcP,sPcP,pPKP,sPKP,pPKIKP,sPKIKP,pS,sS,pSdiff,sSdiff,pScP,sScP,pSKP,sSKP,sSKS,pSKS,pSKIKP,sSKIKP,pSKIKS,sSKIKS'.split(',')
#     CV_PHS_R_1 = CV_PHS_Z_1

#     # surface reflections [1 bounce]
#     CV_PHS_T_2 = CV_PHS_T_1 + 'SS,ScSScS,sSS,sScSScS'.split(',')
#     CV_PHS_Z_2 = CV_PHS_Z_1 + 'PP,SP,SS,PcPPKP,PKPPKP,SKPPKP,SKPSKP,PKIKPPKIKP,pPP,sPP,pSP,sSP,pSS,sSS,pPcPPKP,sPcPPKP,pPKPPKP,sPKPPKP,pSKPPKP,sSKPPKP,pSKPSKP,sSKPSKP,pPKIKPPKIKP,sPKIKPPKIKP'.split(',')
#     CV_PHS_R_2 = CV_PHS_Z_2

#     # surface reflections [2 bounces]
#     CV_PHS_T_3 = CV_PHS_T_2 + 'SSS,sSSS,ScSScSScS,sScSScSScS'.split(',')
#     CV_PHS_Z_3 = CV_PHS_Z_2 + 'PPP,SPP,SSP,SSS,PcPPKPPKP,pPPP,sPPP,pSPP,sSPP,pSSP,sSSP,pSSS,sSSS,pPcPPKPPKP,sPcPPKPPKP'.split(',')
#     CV_PHS_R_3 = CV_PHS_Z_3

#     # surface reflections [3 bounces]
#     CV_PHS_T_4 = CV_PHS_T_3 + 'SSSS,sSSSS'.split(',')
#     CV_PHS_Z_4 = CV_PHS_Z_3 + 'PPPP,pPPPP,sPPPP,SSSS,pSSSS,sSSSS'.split(',')
#     CV_PHS_R_4 = CV_PHS_Z_4


#     # pure P
#     CV_PHS_P_1 = 'p,P,Pdiff,PcP,PKP,PKIKP,pP,sP,pPdiff,pPcP,pPKP,pPKIKP'.split(',')
#     CV_PHS_P_2 = CV_PHS_P_1 + 'PP,PcPPKP,PKPPKP,PKIKPPKIKP,pPP,pPcPPKP,pPKPPKP,pPKIKPPKIKP'.split(',')
#     CV_PHS_P_3 = CV_PHS_P_2 + 'PPP,PcPPKPPKP,PKPPKPPKP,pPPP,pPcPPKPPKP,pPKPPKPPKP'.split(',')
#     CV_PHS_P_4 = CV_PHS_P_3 + 'PPPP,pPPPP'.split(',')

#     CV_PHASES = {'T': CV_PHS_T, 'Z': CV_PHS_Z, 'R': CV_PHS_R}
#     CV_PHASES_1 = {'T': CV_PHS_T_1, 'Z': CV_PHS_Z_1, 'R': CV_PHS_R_1}
#     CV_PHASES_2 = {'T': CV_PHS_T_2, 'Z': CV_PHS_Z_2, 'R': CV_PHS_R_2}
#     CV_PHASES_3 = {'T': CV_PHS_T_3, 'Z': CV_PHS_Z_3, 'R': CV_PHS_R_3}
#     CV_PHASES_4 = {'T': CV_PHS_T_4, 'Z': CV_PHS_Z_4, 'R': CV_PHS_R_4}

#     output_dict = {'Master': CV_PHASES, 'Direct': CV_PHASES_1, '1_refl': CV_PHASES_2, '2_refl': CV_PHASES_3, '3_refl': CV_PHASES_4 }

#     return output_dict


# def main():

#     output_dict = phases()
#     print(output_dict)

# if __name__ == '__main__':
#     main()

