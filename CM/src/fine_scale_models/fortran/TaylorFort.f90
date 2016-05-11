real*8 function norm(tensor)
   real*8 tensor(6)
   real*8 acc
   acc = sqrt(tensor(1)*tensor(1)   &
         + tensor(2)*tensor(2)   &
         + tensor(3)*tensor(3)   &
         + 2.0*tensor(4)*tensor(4)   &
         + 2.0*tensor(5)*tensor(5)   &
         + 2.0*tensor(6)*tensor(6))
   norm = acc
   return 
end


subroutine TaylorF(arr_size, tau_prime, val, m_m, m_D_0, m_g)
   integer arr_size
   real*8 tau_prime(arr_size)
   real*8 val(arr_size)
   real*8 m_m, m_D_0, m_g
   real*8 norm_tau_prime
   real*8 norm

   !write(*,*) 'm_m = ',m_m,' m_D_0 = ', m_D_0,' m_g = ',m_g
   if ( m_m == 1. ) then
      val = m_D_0 * tau_prime / m_g
   else

      norm_tau_prime = norm(tau_prime)
!      write(*,*) 'Norm = ', norm_tau_prime
      !norm_tau_prime = sqrt(tau_prime(1)**2 + tau_prime(2)**2 + tau_prime(3)**2)

      if (norm_tau_prime > 0.0) then
         val = m_D_0 * tau_prime * ( (norm_tau_prime/m_g)**(1.0/m_m) / norm_tau_prime )
      else
         val = 0.0
      end if
   end if
   return
end

subroutine PrintTensor2 (values)
   real*8 values(6)
   write(*,*) (values(i), i=1,6)
   return
end

