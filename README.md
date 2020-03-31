# covid19_nbody
This is an attempt to simulate epidemic spread inspired from 3blue1brown channel.
Here is how the simulation looks:
![Alt Text](./covid19.gif)


# TODO list
* Implement infection radius
* Implement a repulsive force for social distancing
* Implement multiple colonies
* Include demographic information
* Humans dont undergo elastic collisions. The particles should slow down when coming close to other particles but travel with some velocity otherwise.
* Implement quarantine
* Implement a scenario where everyday the severity of disease increases or decreases with some probability.  The probability of death increases with severity
* Include a probability of person getting infected when in contact with infected person as an approximation for better hygiene habits
* Social distancing and hygiene on individual basis. Does social distance and hygiene prove advantageous on individual basis?
* Implement quarantine only if a person is detected positive. The probability of a person being detected depends upon symptoms and symptoms are related to severity. 
* Relate the severity of disease with demographic information. 
