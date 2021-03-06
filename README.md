![balistica](https://user-images.githubusercontent.com/69157203/91104006-d8528600-e63a-11ea-87f4-b0531ee0dc53.png)


# SATELITE SENTINEL 1A/B

El codigo implementado arroja esto en consola.

![image](https://user-images.githubusercontent.com/69157203/92353457-c222dc00-f0b6-11ea-82e4-6e4c39f0815d.png)


Se realiza una prediccion del desplazamineto de un satelite por medio de una ecuacion diferencila de movimiento. 

Se grafican sus posiciones x, y, z por medio de un archivo de texto .EOF 

![SATELITE1](https://user-images.githubusercontent.com/69157203/92353329-69ebda00-f0b6-11ea-9545-8b919accbd69.png)


Se utilizan los valores iniciales obtenidos anterirmente para comparar la integracion numerica bajo dos metodos distintos. ODEINT y EULERINT con una subdivision. 

![SATELITE2](https://user-images.githubusercontent.com/69157203/92353332-6c4e3400-f0b6-11ea-9a0f-7cbb02b8ca02.png)

esto muestra que la solucion comienza a diverger a medida que va aumentando el tiempo, por otro lado, la solucion es menos precisa en eulerint que en odeint y su tiempo de implementacion es mayor para eulerint que para odeint. seindo casi el doble de tiempo para euler int que para odeint.

Se analiza la convergencia de los resultados (Se comentan en los codigos por su nivel de exigencia computacional en cada corrida) y se observa que la diferencia es sinusoidal, lo que tiene señala que para N pequeños, se puede encontrar con alguna solucion particular que este en la vecindad de la solucion y que tenga menos de 1% de error, pero en terminos de precision de la solucion, no indican mucho, ya que difiere en las demás zonas. Para N grandes, se ve que se aplana la amplitud de la sinusioide y se ibserva que tiende a aplanar la solucion. De ser necesario el uso de eulerint para la solucion, seria necesario minimo N = 1000 para poder llegar a una solucion con vecindario razonablemente cercano. 

![Figure 2020-09-07 011516](https://user-images.githubusercontent.com/69157203/92353898-a3711500-f0b7-11ea-9896-3cf240729c77.png)

Se observa que para N=500 se tiene un error cercano a 31% y en su ejecucion se demoro cerca de 200 segundos.

Al agregar las correcciones, se observa un anotable mejora de unos metros, aun cuando ya se tenia una buena precision de los resultados obtenidos en el ciclo asignado.

Se presenta un grafico tridimensional de la orbita del satelite, comprobadose su cercania con la realidad y mejorando cerca de 8 metros en precision, mostrandon una diferencia de 359 m con el modelo.

![SATELITE4](https://user-images.githubusercontent.com/69157203/92354242-580b3680-f0b8-11ea-8316-afec5eb3183e.png)


![SATELITE5](https://user-images.githubusercontent.com/69157203/92354238-56417300-f0b8-11ea-9c25-6297a8273ccb.png)

Demorando un total de 2.5 [s] en correr el codigo.

# Entrega Final
-----------------------------------------------------------------------------------------------------------------------------------------
Se comentan los resultados de la convergencia con la implementacion de los polinomios de Lagendre en un Modelo Geopotencial de la tierra.
-----------------------------------------------------------------------------------------------------------------------------------------

Se implemento el modelo Geopotencial de la tierra, generando el potencial de la tierra. Esto debiese generar una convergencia a la orbita estacionaria del satelite, pero en vez de eso, no se observo cambio. Se generó una forma de incorporar más polinomios de Legendre y coeficientes, usando las funciones de numpy para ello a medida que se agregaban polinomios, no convergia, sino que aumentaba el tiempo de iteracion.

Si bien, el modelo genera una forma similar al del satelite original, no se observa la convergencia con el modelo por si solo, por lo que se concluye que no es suficiente en este caso el modelo geopotencial. 

![SATELITE](https://user-images.githubusercontent.com/69157203/92985428-cd944f80-f488-11ea-926b-8b69a637e627.png)


Se considera la asociacion de las condiciones iniciales, las cuales en un estado estacionario, tiene roce y velocidades asociadas a su inercia y mása,las cuales dificultan su desplazamiento y pudiesen interferir en el modelo. se observa que la convergencia a medida que avanza el tiempo de orbita diverge, asociado a los graficos adjuntos. 

![SATELITE2](https://user-images.githubusercontent.com/69157203/92985544-041e9a00-f48a-11ea-8d7c-c8f9b986c0cf.png)

es probable que la implementacion no haya sido la adecuada del modelo geopotencial porque no se vio cambios significativos, pero no se pudo compensar las diferencias de forma satisfactoria a menos de cambair las condiciones iniciales.


