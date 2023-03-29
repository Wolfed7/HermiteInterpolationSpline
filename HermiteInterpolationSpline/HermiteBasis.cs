namespace HermiteInterpolationSpline
{
   internal class HermiteBasis
   {
      public static double Psi0(double Xi) 
         => 1 - 3 * Xi * Xi + 2 * Xi * Xi * Xi;
      public static double Psi1(double Xi, double h) 
         => h * (Xi - 2 * Xi * Xi + Xi * Xi * Xi);
      public static double Psi2(double Xi) 
         => 3 * Xi * Xi - 2 * Xi * Xi * Xi;
      public static double Psi3(double Xi, double h) 
         => h * (-Xi * Xi + Xi * Xi * Xi);
   }
}
